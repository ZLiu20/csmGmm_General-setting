# For Supp Fig 4B

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/SuppFig9/SFig5B_sim_csmGmm.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig4/Fig4B_sim.R")

# # load libraries
# library(mvtnorm)
# library(data.table)
# library(bindata)
# library(dplyr)
# library(magrittr)
# library(devtools)
# library(ks)
# library(csmGmm)

# load libraries
library(dplyr)
library(mvtnorm)
library(data.table)
library(bindata)
library(magrittr)
library(devtools)
library(rje)
library(ks)
library(csmGmm)
library(here)
library(locfdr)



# define the  function

################################################################################
# define the  function

################################################################################
find_max_means_R1 <- function(muInfo) {
  
  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  S1 <- c(1,2,4,8)+1
  # Com_Set <- compute_sets(K,t)
  # S1 <- Com_Set[1]
  # just keep finding the max
  maxMeans <- rep(0, K) 
  for (element_it in S1) {
    tempMat <- cbind(muInfo[[element_it]], maxMeans)
    maxMeans <- apply(tempMat, 1, max)
  }
  # return K*1 vector
  return(maxMeans)
}


define_H_space_R1 <- function(K,t) {
  H_space <- expand.grid( rep(list(c(-1,0,1)), K) )
  s <- rep(0, nrow(H_space))
  for (col_it in 1:ncol(H_space)) {
    s <- s + abs(H_space[, col_it])
  }
  H_space <- H_space %>% mutate(s = s) %>%
    arrange(s) %>%
    select(-s) %>%
    as.matrix(.)
  
  H_annot <- rep(0, nrow(H_space))
  # H_annot[which(apply(abs(H_space), 1, sum) == K)] <- 1
  H_annot[which(t+1<=apply(abs(H_space), 1, sum) & apply(abs(H_space), 1, sum) <= K)] <- 1
  
  
  return( list(H_space=H_space, H_annot=H_annot) )
}


################################################################################
emp_bayes_framework_R1 <- function(t_value, summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE,
                                   ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
  # Step 1 
  H_space_output <- define_H_space_R1(K = ncol(summary_tab), t=t_value)
  H_space <- as.matrix(H_space_output$H_space)
  H_annot <- H_space_output$H_annot
  if (sameDirAlt) {
    #null_rows <- which(apply(H_space, 1, sum) != ncol(summary_tab) & apply(H_space, 1, sum) != -ncol(summary_tab))
    null_rows <- which(apply(abs(H_space), 1, sum) <= t_value |
                         (apply(abs(H_space), 1, sum) > t_value &
                            apply(abs(H_space), 1, sum) <= ncol(summary_tab) &
                            apply(abs(H_space), 1, sum) != abs(apply(sign(H_space), 1, sum)))) #revised on 15/11/2025
  } else { 
    null_rows <- which(H_annot == 0)
  }
  
  ###################################
  # Step 2
  # Find the distribution f and pi(H)
  locfdr_output <- locfdr_estim(summary_tab=summary_tab, kernel=kernel, joint=joint, ind=ind, dfFit = dfFit, nulltype = nulltype)
  
  # can be negative - see first real data analysis
  pi0_estim <- locfdr_output$pi0_estim
  if (length(which(pi0_estim >= 1 | pi0_estim <= 0)) > 0) {
    pi0_estim[which(pi0_estim >= 1 | pi0_estim <= 0)] <- 0.99
  }
  marg_pmf_tab <- locfdr_output$marg_pmf_tab
  # was an error with locfdr
  if (is.null(marg_pmf_tab)) {return(-1)} 
  if (ncol(marg_pmf_tab) < 2*ncol(summary_tab)) {
    return(-1)
  }
  
  # calculate conditional distributions
  H_dist <- calc_H_dist_indep(pi0_estim=pi0_estim, H_space=H_space)
  cond_pmfs <- calc_cond_pmfs(summary_tab=summary_tab, marg_pmf_tab=marg_pmf_tab, pi0_estim=pi0_estim)
  null_pmf_tab <- cond_pmfs$null_pmf_tab
  neg_pmf_tab <- cond_pmfs$neg_pmf_tab
  pos_pmf_tab <- cond_pmfs$pos_pmf_tab
  
  # check
  if (is.na(sum(neg_pmf_tab)) | is.na(sum(neg_pmf_tab))) {return(-1)}
  
  # run EM
  binned_tab <- cond_pmfs$binned_tab
  EMoutput <- run_EM_pmf(binned_dat=binned_tab, H_dist=H_dist, H_space=H_space,
                         null_pmf_tab=null_pmf_tab, neg_pmf_tab=neg_pmf_tab,
                         pos_pmf_tab=pos_pmf_tab, pi0_estim=pi0_estim, epsilon=Hdist_epsilon)
  
  ###########################
  # Step 3: Calculate the local Bayes FDR
  loc_fdr_num <- apply(EMoutput$jointMat[, null_rows], 1, sum)
  lfdrVec <- loc_fdr_num / EMoutput$probZ
  
  return(list(marg_pmf_tab = marg_pmf_tab, null_pmf_tab = null_pmf_tab, pi0_estim = pi0_estim, Hdist_final = EMoutput$H_dist,
              neg_pmf_tab = neg_pmf_tab, pos_pmf_tab = pos_pmf_tab, lfdrVec = lfdrVec))
}


################################################################################
symm_fit_ind_EM_R1 <- function(t_value, testStats, initMuList, initPiList, sameDirAlt=FALSE, eps = 10^(-5), checkpoint=TRUE) {
  
  # number of composite null hypotheses
  J <- nrow(testStats)
  # number of dimensions
  K <- ncol(testStats)
  B <- 2^K - 1
  # number of hl configurations - 1
  L <- 3^K - 1
  #number of the non_zero values
  t <- t_value
  # make all configurations
  Hmat <- expand.grid(rep(list(-1:1), K))
  # attach the bl
  blVec <- rep(0, nrow(Hmat))
  slVec <- rep(0, nrow(Hmat))
  for (k_it in K:1) {
    blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
    slVec <- slVec + abs(Hmat[, k_it]) #number of non_zeros
  }
  # symmetric alternative
  sum_Hmat <- apply(Hmat[, 1:K], 1, sum)
  symAltVec <- ifelse((t+1 <= sum_Hmat & sum_Hmat <= K) |  (sum_Hmat >= -K & sum_Hmat <= -t-1), 1, 0)
  #symAltVec <- ifelse(apply(Hmat[, 1:K], 1, sum) == K | apply(Hmat[, 1:K], 1, sum) == -K, 1, 0)
  
  # sort Hmat
  Hmat <- Hmat %>% dplyr::mutate(bl = blVec) %>%
    dplyr::mutate(sl = slVec) %>%
    dplyr::mutate(symAlt = symAltVec) %>%
    dplyr::arrange(.data$bl, .data$Var1, .data$Var2) %>%
    dplyr::mutate(l = 0:(nrow(.) - 1)) %>%
    dplyr::relocate(.data$l, .before = .data$symAlt)
  
  # initialize
  # the muInfo and piInfo are lists that hold the information in a more compact manner.
  # the allMu and allPi are matrices that repeat the information so the calculations can be performed faster.
  muInfo <- initMuList
  piInfo <- initPiList
  oldParams <- c(unlist(piInfo), unlist(muInfo))
  MbVec <- sapply(piInfo, FUN=length)
  
  # run until convergence
  diffParams <- 10
  iter <- 0
  while (diffParams > eps) {
    
    #####################################################
    # first, update allPi and allMu with the new parameter values
    
    # allPi holds the probabilities of each configuration (c 1), the bl of that
    # configuration (c 2), the sl of that configuration (c3), the m of that configuration (c 4),
    # the l of that configuration (c 5), and whether it's part of the symmetric alternative (c 6).
    
    # number of rows in piMat is L if Mb = 1 for all b.
    # number of rows is \sum(l = 0 to L-1) {Mbl}
    allPi <- c(piInfo[[1]], 0, 0, 1, 0, 0)
    
    # allMu holds the mean vectors for each configuration in each row, number of columns is K
    allMu <- rep(0, K)
    
    # loop through possible values of bl
    for (b_it in 1:B) {
      
      # Hmat rows with this bl
      tempH <- Hmat %>% dplyr::filter(.data$bl == b_it)
      
      # loop through possible m for this value of bl
      for (m_it in 1:MbVec[b_it + 1]) {
        allPi <- rbind(allPi, cbind(rep(piInfo[[b_it + 1]][m_it] / (2^tempH$sl[1]), nrow(tempH)), tempH$bl,
                                    tempH$sl, rep(m_it, nrow(tempH)), tempH$l, tempH$symAlt))
        for (h_it in 1:nrow(tempH)) {
          allMu <- cbind(allMu, unlist(tempH %>% dplyr::select(-.data$bl, -.data$sl, -.data$l, -.data$symAlt) %>%
                                         dplyr::slice(h_it)) * muInfo[[b_it + 1]][, m_it])
        }
      } # done looping through different m
    } # dont looping through bl
    colnames(allPi) <- c("Prob", "bl", "sl", "m", "l", "symAlt")
    
    ##############################################################################################
    # this is the E step where we calculate Pr(Z|c_l,m) for all c_l,m.
    # each l,m is one column.
    # only independence case for now
    if (ncol(testStats) == 2) {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_2d, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    } else if (ncol(testStats) == 3) {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_3d, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    } else {
      conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_ind_multiple, Zmat = testStats) %>%
        sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    }
    probZ <- apply(conditionalMat, 1, sum)
    AikMat <- conditionalMat %>% sweep(., MARGIN=1, STATS=probZ, FUN="/")
    Aik_alln <- apply(AikMat, 2, sum) / J
    
    ###############################################################################################
    # this is the M step for probabilities of hypothesis space
    for (b_it in 0:B) {
      for (m_it in 1:MbVec[b_it + 1]) {
        tempIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        piInfo[[b_it + 1]][m_it] <- sum(Aik_alln[tempIdx])
      }
    }
    
    # M step for the means
    # loop through values of bl
    # do the alternative last, enforce that it must be larger in magnitude than the nulls
    for (b_it in 1:B) {
      
      tempHmat <- Hmat %>% dplyr::filter(.data$bl == b_it)
      # loop through m
      for (m_it in 1:MbVec[b_it + 1]) {
        tempMuSum <- rep(0, nrow(allMu))
        tempDenom <- rep(0, nrow(allMu))
        
        # these are the classes that contribute to \mu_bl,m
        AikIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        for (idx_it in 1:length(AikIdx)) {
          tempAik <- AikIdx[idx_it]
          tempHvec <- tempHmat %>% dplyr::select(-.data$bl, -.data$sl, -.data$l, -.data$symAlt) %>%
            dplyr::slice(idx_it) %>% unlist(.)
          
          tempMuSum <- tempMuSum + colSums(AikMat[, tempAik] * sweep(x = testStats, MARGIN = 2,
                                                                     STATS = tempHvec, FUN="*"))
          tempDenom <- tempDenom + rep(J * Aik_alln[tempAik], length(tempDenom)) * abs(tempHvec)
        } # done looping for one l, m
        whichZero <- which(tempDenom == 0)
        tempDenom[whichZero] <- 1
        muInfo[[b_it + 1]][, m_it] <- tempMuSum / tempDenom
        
        # make sure mean constraint is satisfied
        S2 <- c(3,5,6,7,9,10,11,12,13,14,15)
        # Com_Set <- compute_sets(K,t)
        # S2 <- Com_Set[2]
        
        if (b_it %in% S2) {
          maxMeans <- find_max_means_R1(muInfo)
          whichSmaller <- which(muInfo[[b_it+1]][, m_it] < maxMeans)
          if (length(whichSmaller) > 0) {
            muInfo[[b_it+1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
            #if (b_it == B) {
            
            #maxMeans <- find_max_means_R1(muInfo)
            #whichSmaller <- which(muInfo[[b_it + 1]][, m_it] < maxMeans)
            #if (length(whichSmaller) > 0) {
            #  muInfo[[b_it + 1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
            
          }
        } # done with mean constraint
        
      } # done looping through m
      
    } # done updating means
    
    ###############################################################################################
    # find difference
    allParams <- c(unlist(piInfo), unlist(muInfo))
    diffParams <- sum((allParams - oldParams)^2)
    
    # update
    oldParams <- allParams
    iter <- iter + 1
    if (checkpoint) {
      cat(iter, " - ", diffParams, "\n", allParams, "\n")
    }
  }
  
  # calculate local fdrs
  if (sameDirAlt) {
    nullCols <- which(allPi[, 6] == 0)
  } else {
    nullCols <- which(allPi[, 3] <= t)  # revised on 2025/09/28
    #nullCols <- which(allPi[, 3] < K)
  }
  probNull <- apply(conditionalMat[, nullCols], 1, sum)
  lfdrResults <- probNull / probZ
  
  return(list(piInfo = piInfo, muInfo = muInfo, iter = iter,
              lfdrResults = lfdrResults))
}



################################################################################
check_incongruous_R1 <- function(zMatrix, lfdrVec, t_value) {
  K <- ncol(zMatrix)
  
  # Step 1: Define H space
  H_output <- define_H_space_R1(K = K, t = t_value)
  H_space <- H_output$H_space
  H_annot <- H_output$H_annot
  
  # Step 2: Assign configuration to each SNP
  sign_mat <- sign(zMatrix)
  config_idx <- apply(sign_mat, 1, function(s) {
    match_idx <- which(apply(H_space, 1, function(h) all(h == s)))
    if (length(match_idx) == 0) return(NA)  # Z=0 
    return(match_idx[1])
  })
  
  # Remove SNPs with undefined config (e.g., all Z=0)
  valid <- !is.na(config_idx)
  if (!any(valid)) return(c())
  
  zMatrix <- zMatrix[valid, , drop = FALSE]
  lfdrVec <- lfdrVec[valid]
  config_idx <- config_idx[valid]
  is_alt <- H_annot[config_idx]
  
  # Step 3: Get alternative and null sets
  alt_idx <- which(is_alt == 1 & lfdrVec < 0.99)
  null_idx <- which(is_alt == 0 & lfdrVec < 0.99)
  
  if (length(alt_idx) == 0 || length(null_idx) == 0) return(c())
  
  z_alt <- abs(zMatrix[alt_idx, , drop = FALSE])
  lfdr_alt <- lfdrVec[alt_idx]
  z_null <- abs(zMatrix[null_idx, , drop = FALSE])
  lfdr_null <- lfdrVec[null_idx]
  
  # Step 4: Check inconsistency
  bad_alt <- c()
  for (i in seq_len(nrow(z_alt))) {
    # Find null points with smaller magnitude in ALL dimensions
    dominated <- apply(z_null, 1, function(zn) all(zn < z_alt[i, ]))
    if (any(dominated)) {
      # If any dominated null point has lfdr <= current alt point ? incongruous
      if (any(lfdr_null[dominated] <= lfdr_alt[i])) {
        bad_alt <- c(bad_alt, which(valid)[alt_idx[i]])  # map back to original index
      }
    }
  }
  
  return(bad_alt)
}


################################################################################
check_incongruous <- function(zMatrix, lfdrVec) {
  
  # remove lfdr = 1
  lessThanOne <- which(lfdrVec < 0.99)
  if (length(lessThanOne) <= 1) {return(c())}
  zMatrix <- zMatrix[lessThanOne, ]
  lfdrVec <- lfdrVec[lessThanOne]
  
  # do it in K^2 quadrants
  K <- ncol(zMatrix)
  quadrants <- expand.grid(rep(list(c(-1, 1)), K))
  
  badIdx <- c()
  for (quad_it in 1:nrow(quadrants)) {
    # separate into quadrants
    idxVec <- 1:nrow(zMatrix)
    tempStats <- zMatrix
    tempLfdr <- lfdrVec
    for (k_it in 1:K) {
      if (class(tempStats)[1] == "numeric") {break}
      if (quadrants[quad_it, k_it] == -1) {
        toKeep <- which(tempStats[, k_it] < 0 )
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      } else {
        toKeep <- which(tempStats[, k_it] > 0 )
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      }
      tempStats <- tempStats[toKeep, ]
    } # end finding quadrant
    if (length(idxVec) <= 1) {next}
    # take absolute value
    tempStats <- abs(tempStats)
    
    # order by lfdr
    tempDat <- tempStats %>% as.data.frame(.data) %>%
      dplyr::mutate(lfdr = tempLfdr) %>%
      dplyr::mutate(idx = idxVec)
    
    # check for incongruous
    if (K == 2) {
      colnames(tempDat)[1:2] <- c("Z1", "Z2")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2))
      incongruousVec <- sapply(1:nrow(tempDat),FUN = find_2d,  allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2)))
    } else if (K == 3) {
      colnames(tempDat)[1:3] <- c("Z1", "Z2", "Z3")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2), dplyr::desc(.data$Z3))
      incongruousVec <- sapply(1:nrow(tempDat), FUN = find_3d, allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2, .data$Z3)))
    } else {
      stop("only support for 2-3 dimensions right now")
    }
    
    # get the bad indices
    badIdx <- c(badIdx, tempDat$idx[which(incongruousVec > 0)])
  }
  
  return(badIdx)
}


#' Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
#' test statistics with lower magnitudes). For K=2 case.
#'
#' @param x Scalar, which row of allTestStats to check.
#' @param allTestStats J*K vector of all test statistics.
#'
#' @return A scalar denoting the number of sets with lower lfdr and test statistics of lower magnitude. 0 means congruous result.
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4))
#' find_2d(x = 5, allTestStats = zMatrix)
#'
find_2d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2]))
}

#' Tells if row x if allTestStats is an incongruous result (has a higher lfdr than a set of
#' test statistics with lower magnitudes). For K=3 case.
#'
#' @param x Scalar, which row of allTestStats to check.
#' @param allTestStats J*K vector of all test statistics.
#'
#' @return A scalar denoting the number of sets with lower lfdr and test statistics of lower magnitude. 0 means congruous result.
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4),  rnorm(10^4))
#' find_3d(x = 5, allTestStats = zMatrix)
#'
find_3d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2] &
                 allTestStats[1:x, 3] < allTestStats[x, 3]))
}

################################################################################



# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Fig4", "output")
outName <- paste0(outputDir, "/Fig4B_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- FALSE
# the name will be [testStatsName]_[betaStart]_S[Snum]_aID[aID].txt
testStatsName <- here::here(outputDir, "allZ")
betaName <- here::here(outputDir, "allBeta")

# parameters
doHDMT <- FALSE
doDACT <- FALSE
doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE
doNew <- TRUE
qvalue <- 0.1
nSNPs <- 10^5 
setSize <- 1000
n <- 1000
nDims <- 4
t <- 2
nSets <- nSNPs / setSize
nSims <- 5
margprob <- rep(0.3, setSize)
simsPerEffSize <- 20
effSizeMult <- ceiling(aID / simsPerEffSize)
# betaStart <- 0.24
# betaMin <- rep(betaStart + 0.01 * effSizeMult, nDims)
# betaMax <- rep(betaStart + 0.1 + 0.01 * effSizeMult, nDims)
## revised on 2025/12/09
# betaMin <- rep(0.14 + 0.01 * effSizeMult, nDims)
# betaMax <- rep(0.24 + 0.1 + 0.01 * effSizeMult, nDims)
# betaMin <- rep(0.10 + 0.01 * effSizeMult, nDims)
# betaMax <- rep(0.20 + 0.1 + 0.01 * effSizeMult, nDims)
# betaMin <- rep(0.05 + 0.01 * effSizeMult, nDims)
# betaMax <- rep(0.20 + 0.1 + 0.01 * effSizeMult, nDims)
betaMin <- rep(0.20 + 0.01 * effSizeMult, nDims)
betaMax <- rep(0.20 + 0.1 + 0.01 * effSizeMult, nDims)
beta0 <- -1

# determines how many signals there are
sProp <- c(0.9579, 0.04, 0.0012, 0.0006, 0.0003)
hMat <- expand.grid(rep(list(c(-1, 0, 1)), nDims)) %>%
  as.matrix(.) %>%
  cbind(., rowSums(abs(.))) %>%
  as.data.frame(.) %>%
  set_colnames(c(paste0("Var", 1:(ncol(.)-1)), "s")) %>%
  arrange(s)
number <- c()
for (s_it in 0:max(hMat$s)) {
  numRows <- length(which(hMat$s == s_it))
  number <- c(number, rep(sProp[s_it + 1] * nSNPs / numRows, numRows))
}
hMat <- hMat %>% mutate(number = number) %>%
  mutate(number = round(number))

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  minEff1=betaMin[1],
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)
# each loop is one simulation iteration
for (sim_it in 1:nSims) {
  
  # set the seed 
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it
  
  # load or save data
  if (loadData) {
    allZ <- fread(paste0(testStatsName, "_", betaStart, "_S", Snum, "_aID", aID, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_", betaStart, "_S", Snum, "_aID", aID, ".txt"), data.table=F)
  } else {
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    
    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)
    
    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets)  {
      # save test statistics      
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)
      
      # generate coefficient matrix
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                          betaMin=betaMin, betaMax=betaMax)
      
      # loop through each dimension
      for (dimension_it in 1:nDims) {
        tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
        tempCoef <- coefMat[, dimension_it]
        
        # outcome
        tempEta <- sweep(tempG, MARGIN=2, STATS=tempCoef, FUN="*") + matrix(data=beta0, nrow=nrow(tempG), ncol=ncol(tempG))
        tempMu <- rje::expit(as.numeric(tempEta))
        tempY <- rbinom(n=length(tempMu), size=1, prob=tempMu)
        # put back in matrix form
        yMat <- matrix(data=tempY, nrow=nrow(tempEta), ncol=ncol(tempEta), byrow=FALSE)
        
        # loop through each SNP
        for (test_it in 1:ncol(yMat)) {
          tempMod <- glm(yMat[, test_it] ~ tempG[, test_it], family=binomial)
          statsMat[test_it, dimension_it] <- summary(tempMod)$coefficients[2, 3]
        }
        cat("simulated dimension ", dimension_it, "\n")
      }
      
      # record data
      startIdx <- (set_it - 1) * setSize + 1
      endIdx <- set_it * setSize
      allZ[startIdx:endIdx, ] <- statsMat
      allBeta[startIdx:endIdx, ] <- coefMat
      
      # checkpoint 
      if(set_it%%1 == 0) {cat(set_it)}
    } # done generating data
    
    # save it
    if (saveData) { 
      write.table(allZ, paste0(testStatsName, "_", betaStart, "_S", Snum, "_aID", aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_", betaStart, "_S", Snum, "_aID", aID, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    } 
  }
  
  # number of signals and causal SNPs
  #causalVec <- as.numeric(apply(allBeta, 1, prod) != 0) 
  ## revised on 30/11/2025
  causalVec <- as.numeric(rowSums(allBeta != 0) >= t+1 & rowSums(allBeta != 0) <= nDims)
  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$pi0aTrue[sim_it] <- length(which(allBeta[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allBeta[, 2] != 0))
  
  # adjustment to not get p-values of 0 needed for DACT and HDMT 
  for (col_it in 1:ncol(allZ)) {
    tooBig <- which(allZ[, col_it] > 8.1)
    tooSmall <- which(allZ[, col_it] < -8.1)
    if (length(tooBig) > 0) {
      allZ[tooBig, col_it] <- 8.1
    }
    if (length(tooSmall) > 0) {
      allZ[tooSmall, col_it] <- -8.1
    }
  }
  # p-value matrix
  allP <- 1- pchisq(as.matrix(allZ)^2, df=1)
  
  # hold the results
  totOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], origIdx = 1:nrow(allP), causal=causalVec) %>%
    mutate(pmax = pmax(X1, X2))
  
  # analyze it 
  # start with HDMT
  nullprop <- tryCatch(null_estimation(allP), error=function(e) e, warning=function(w) w)
  if (class(nullprop)[1] == "list" & doHDMT) {
    hdmtRes <- tryCatch(fdr_est_orig(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,allP),
                        error = function(e) e, warning = function(w) w)
    if (class(hdmtRes)[1] == "numeric") {
      
      totOut <- totOut %>% mutate(hdmtRes = hdmtRes)
      
      # calculate power and fdp for HDMT
      powerRes$fdpHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$causal == 0)) /  length(which(totOut$hdmtRes < qvalue))
      powerRes$powerHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejHDMT[sim_it] <- length(which(totOut$hdmtRes < qvalue)) 
    } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  
  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  if (doDACT) {
    DACTout <- tryCatch(DACT_noEst(p_a = allP[, 1], p_b = allP[, 2], nullEst=nullprop, correction="JC"),
                        error = function(e) e, warning=function(w) w)
    if (class(DACTout)[1] %in% c("simpleError", "simpleWarning")) {
      totOut <-  totOut %>% mutate(DACTp = NA)
    } else {
      DACTdf <- data.frame(DACTp = DACTout$pval, origIdx = 1:length(DACTout$pval)) %>%
        arrange(DACTp) %>%
        mutate(rankedIdxP = 1:nrow(.)) %>%
        mutate(km = 1:nrow(.)/ nrow(.)) %>%
        mutate(RHS = km * qvalue) 
      rejected <- which(DACTdf$DACTp <= DACTdf$RHS)
      if (length(rejected) == 0) {
        maxIdx <- 0
      } else {maxIdx <- max(rejected)}
      DACTdf <- DACTdf %>% mutate(reject = ifelse(rankedIdxP <= maxIdx, 1, 0)) %>%
        arrange(origIdx)
      # append results
      totOut <- totOut %>% mutate(DACTp = DACTdf$DACTp, DACTrej = DACTdf$reject)
      
      # power and FDP for DACT
      powerRes$fdpDACT[sim_it] <- length(which(totOut$DACTrej == 1 & totOut$causal == 0)) / length(which(totOut$DACTrej == 1))
      powerRes$powerDACT[sim_it] <- length(which(totOut$DACTrej == 1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejDACT[sim_it] <- length(which(totOut$DACTrej == 1)) 
    }
  } else {totOut <-  totOut %>% mutate(DACTp = NA, DACTrej = NA)}
  
  # old method - kernel
  if (doKernel) { 
    # oldResKernel <- emp_bayes_framework(summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE, 
    #                                     dfFit = 7, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE,
                                           dfFit = 7, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 0)) / length(which(totOut$kernelAvg < qvalue))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
      # powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$inconKernel[sim_it] <- length(check_incongruous_R1(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec,t_value = t))
      powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue)) 
    } else {
      totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)
    }
  } else {totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)}
  
  # old method - 7 df
  if (do7df) { 
    # oldRes7df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
    #                                  dfFit = 7, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                        dfFit = 7, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    if (class(oldRes7df)[1] == "list") {
      totOut <- totOut %>% mutate(df7Lfdr = oldRes7df$lfdrVec) %>%
        arrange(df7Lfdr) %>%
        mutate(df7Avg = cummean(df7Lfdr)) %>% 
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdp7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$causal == 0)) /  length(which(totOut$df7Avg < qvalue))
      powerRes$power7df[sim_it] <- length(which(totOut$df7Avg < qvalue & totOut$causal == 1)) / sum(causalVec)
      # powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
      powerRes$incon7df[sim_it] <- length(check_incongruous_R1(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec,t_value = t))
      powerRes$nRej7df[sim_it] <- length(which(totOut$df7Avg < qvalue))
    } else {
      totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)
    } 
  } else {totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)}
  
  # old method - 50 df
  if (do50df) {
    # oldRes50df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
    #                                   dfFit = 50, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                         dfFit = 50, Hdist_epsilon=10^(-1), checkpoint=TRUE)
    if (class(oldRes50df)[1] == "list") {
      totOut <- totOut %>% mutate(df50Lfdr = oldRes50df$lfdrVec) %>%
        arrange(df50Lfdr) %>%
        mutate(df50Avg = cummean(df50Lfdr)) %>% 
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdp50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$causal == 0)) /  length(which(totOut$df50Avg < qvalue))
      powerRes$power50df[sim_it] <- length(which(totOut$df50Avg < qvalue & totOut$causal == 1)) / sum(causalVec)
      # powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
      powerRes$incon50df[sim_it] <- length(check_incongruous_R1(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec,t_value = t))
      powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < qvalue)) 
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}
  
  # new method
  if (doNew) {
    # timing
    startTimeNew <- Sys.time()
    
    initPiList <- list(c(0.82))
    for (i in 2:(2^nDims)) {initPiList[[i]] <- 0.18 / (2^nDims - 1)}
    # the symm_fit_ind.R code will add the appropriate 0s to initMuList
    initMuList <- list(matrix(data=rep(0, nDims), nrow=nDims, ncol=1))
    for (i in 2:(2^nDims)) {
      initMuList[[i]] <- matrix(data=rep(3, nDims), nrow=nDims, ncol=1)
    }
    #newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-1))
    newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))
    
    # timing
    endTimeNew <- Sys.time()
    totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
      arrange(newLfdr) %>%
      mutate(newAvg = cummean(newLfdr)) %>% 
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 0)) /  length(which(totOut$newAvg < qvalue))
    powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
    # powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = newRes$lfdrResults))
    powerRes$inconNew[sim_it] <- length(check_incongruous_R1(zMatrix = allZ, lfdrVec = newRes$lfdrResults,t_value = t))
    powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < qvalue)) 
    powerRes$timeNew[sim_it] <- endTimeNew - startTimeNew 
    
  }
  cat('\n Done with ', sim_it, '\n')
}

write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')

