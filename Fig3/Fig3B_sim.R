# For Figure 3B

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/Fig1B_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig3/Fig3B_sim.R")

# load libraries
library(dplyr)
library(mvtnorm)
library(data.table)
library(bindata)
library(dplyr)
library(magrittr)
library(rje)
library(ks)
library(csmGmm)
library(here)
library(locfdr)


# define the  function
# Function to compute sets N and C
compute_sets <- function(K, t) {
  # Initialize sets N and C
  N <- integer(0)
  C <- integer(0)
  
  # Generate all possible configurations
  for (i in 0:(3^K - 1)) {
    # Convert to base 3 to represent configurations
    base3 <- as.integer(intToBase3(i, K))
    
    # Calculate the number of non-zero elements
    non_zero_count <- sum(abs(base3) == 1)
    
    # Calculate b_l based on the definition
    b_l <- sum(2^(K - seq_len(K)) * (base3 == 1))
    
    # Classify into N or C based on the count of non-zero elements
    if (non_zero_count <= t) {
      N <- unique(c(N, b_l))
    } else {
      C <- unique(c(C, b_l))
    }
  }
  
  # Return sorted unique sets
  return(list(N = sort(unique(N)), C = sort(unique(C))))
}

# Helper function to convert an integer to base 3
intToBase3 <- function(n, length) {
  base3 <- integer(length)
  for (j in seq(length, 1)) {
    base3[j] <- n %% 3
    n <- n %/% 3
  }
  return(base3)
}

##
find_max_means_R1 <- function(muInfo) {
  
  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  S1 <- c(1,2,4)+1
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


emp_bayes_framework_R1 <- function(t_value, summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE,
                                   ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
  # Step 1 
  H_space_output <- define_H_space_R1(K = ncol(summary_tab), t=t_value)
  H_space <- as.matrix(H_space_output$H_space)
  H_annot <- H_space_output$H_annot
  if (sameDirAlt) {
    null_rows <- which(apply(H_space, 1, sum) != ncol(summary_tab) & apply(H_space, 1, sum) != -ncol(summary_tab))
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
        S2 <- c(3,5,6,7)
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
#' check_incongruous.R
#'
#' Check the number of sets of test statistics that have a higher (less significant) lfdr value
#' than other sets with test statistics of uniformly smaller magnitudes.
#'
#' @param zMatrix J*K vector of all test statistics.
#' @param lfdrVec J*1 vector of lfdr values corresponding to each set of test statistics.
#'
#' @return A vector with all the indices of all sets that have a higher lfdr value those a set
#' with smaller test statistic magnitudes.
#' @importFrom dplyr mutate arrange desc select
#' @importFrom rlang .data
#'
#' @export
#' @examples
#' zMatrix <- cbind(rnorm(10^4), rnorm(10^4))
#' lfdrVec <- runif(10^4)
#' check_incongruous(zMatrix = zMatrix, lfdrVec = lfdrVec)
#'
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
outputDir <- here::here("Fig3", "output")
outName <- paste0(outputDir, "/Fig3B_aID", aID, ".txt")

# option to save or load intermediate data to save time,
# set as FALSE for first run and then TRUE thereafter
loadData <- FALSE
saveData <- FALSE
# these names are for if saveData <- TRUE
testStatsName <- here::here(outputDir, "Fig3B_allZ")
betaName <- here::here(outputDir, "Fig3B_allBeta")

# simulation parameters
# doHDMT <- TRUE
# doDACT <- TRUE
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
nDims <- 3
nSets <- nSNPs / setSize
nSims <- 10
t <- 1
margprob <- rep(0.3, setSize)
simsPerEffSize <- 20
effSizeMult <- ceiling(aID / simsPerEffSize)
betaMin <- c(0.38, 0.38, 0.38)
betaMax <- c(0.38, 0.38, 0.38)
# betaMin <- c(0.38, 0.38, 0.38) #result is good, but power is not enough
# betaMax <- c(0.58, 0.58, 0.58)
# betaMin <- c(0.38, 0.38, 0.38)
# betaMax <- c(0.55, 0.55, 0.55)
beta0 <- -1

# determines how many signals there are
pi1 <- 0.01 * effSizeMult
pi11 <- 3 * pi1^2
pi111 <- pi11 / 2
pi00 <- 1 - 3 * pi1 - 3 * pi11 - pi111
sProp <- c(pi00, 3 * pi1, 3 * pi11, pi111)
hMat <- expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)) %>%
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
powerRes <- data.frame(nCausal=rep(NA, nSims), minEff1=pi1,
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA,
                       nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA, nRejNew=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)
for (sim_it in 1:nSims) {
  
  # set the seed 
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it

  # load or save data
  if (loadData) {
    setwd(outputDir)
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else { 
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=nDims)

    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims=nDims)
  
    # generate data in multiple sets - faster than all at once 
    for (set_it in 1:nSets)  {

      # information we need to save
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)

      # generate coefficient matrix    
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                       betaMin=betaMin, betaMax=betaMax) 
  
      ########################################### 
      # three dimensional pleiotropy
      tempG1 <-  sapply(X=margprob, FUN=rbinom, n=n, size=2)
      tempG2 <-  sapply(X=margprob, FUN=rbinom, n=n, size=2)
      tempG3 <- sapply(X=margprob, FUN=rbinom, n=n, size=2) 
      tempAlpha <- coefMat[, 1]
      tempBeta <- coefMat[, 2]
      tempGamma <- coefMat[, 3]
    
      # generate first outcome
      tempEta1 <- sweep(tempG1, MARGIN=2, STATS=tempAlpha, FUN="*") + matrix(data=beta0, nrow=nrow(tempG1), ncol=ncol(tempG1))
      tempMu1 <- rje::expit(as.numeric(tempEta1))
      tempY1 <- rbinom(n=length(tempMu1), size=1, prob=tempMu1)
      yMat1 <- matrix(data=tempY1, nrow=nrow(tempEta1), ncol=ncol(tempEta1), byrow=FALSE)
      # generate second outcome 
      tempEta2 <- sweep(tempG2, MARGIN=2, STATS=tempBeta, FUN="*") + matrix(data=beta0, nrow=nrow(tempG2), ncol=ncol(tempG2))
      tempMu2 <- rje::expit(as.numeric(tempEta2))
      tempY2 <- rbinom(n=length(tempMu2), size=1, prob=tempMu2)
      yMat2 <- matrix(data=tempY2, nrow=nrow(tempEta2), ncol=ncol(tempEta2), byrow=FALSE)
      # generate third outcome
      tempEta3 <- sweep(tempG3, MARGIN=2, STATS=tempGamma, FUN="*") + matrix(data=beta0, nrow=nrow(tempG3), ncol=ncol(tempG3))
      tempMu3 <- rje::expit(as.numeric(tempEta3))
      tempY3 <- rbinom(n=length(tempMu3), size=1, prob=tempMu3)
      yMat3 <- matrix(data=tempY3, nrow=nrow(tempEta3), ncol=ncol(tempEta3), byrow=FALSE)

      # calculate all test statistics 
      for (test_it in 1:ncol(yMat1)) {
        tempMod1 <- glm(yMat1[, test_it] ~ tempG1[, test_it], family=binomial)
        tempMod2 <- glm(yMat2[, test_it] ~ tempG2[, test_it], family=binomial)
        tempMod3 <- glm(yMat3[, test_it] ~ tempG3[, test_it], family=binomial)
        statsMat[test_it, 1] <- summary(tempMod1)$coefficients[2, 3]
        statsMat[test_it, 2] <- summary(tempMod2)$coefficients[2, 3]
        statsMat[test_it, 3] <- summary(tempMod3)$coefficients[2, 3]
      } 
   
      # record
      startIdx <- (set_it - 1) * setSize + 1
      endIdx <- set_it * setSize
      allZ[startIdx:endIdx, ] <- statsMat
      allBeta[startIdx:endIdx, ] <- coefMat
   
      # checkpoint 
      if(set_it%%1 == 0) {cat(set_it)}
    } # done generating data

    # save it
    if (saveData) { 
      setwd(outputDir)
      write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    }
  }

  # number of signals and causal SNPs
  #causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0 & allBeta[, 3] != 0)  # revised on 30/11/2025
  causalVec <- as.numeric(rowSums(allBeta != 0) >= t+1 & rowSums(allBeta != 0) <= nDims)
  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$pi0aTrue[sim_it] <- length(which(allBeta[, 1] != 0))
  powerRes$pi0bTrue[sim_it] <- length(which(allBeta[, 2] != 0))
  powerRes$pi0cTrue[sim_it] <- length(which(allBeta[, 3] != 0))

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
  totOut <- data.frame(X1 = allP[, 1], X2 = allP[, 2], X3 = allP[, 3], origIdx = 1:nrow(allP), causal=causalVec)
 
  # old method - kernel
  if (doKernel) {
    oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE,
                                        dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% 
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 0)) / length(which(totOut$kernelAvg < 0.1))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1)) 
    } else {
      totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)
    }
  } else {totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)}


  # old method - 7 df
  if (do7df) {
    oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldRes7df)[1] == "list") {
      totOut <- totOut %>% mutate(df7Lfdr = oldRes7df$lfdrVec) %>%
        arrange(df7Lfdr) %>%
        mutate(df7Avg = cummean(df7Lfdr)) %>%
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdp7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df7Avg < 0.1))
      powerRes$power7df[sim_it] <- length(which(totOut$df7Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
      powerRes$nRej7df[sim_it] <- length(which(totOut$df7Avg < 0.1))
    } else {
      totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df7Lfdr = NA, df7Avg=NA)}

  # old method - 50 df
  if (do50df) {
    oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                      dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    if (class(oldRes50df)[1] == "list") {
      totOut <- totOut %>% mutate(df50Lfdr = oldRes50df$lfdrVec) %>%
        arrange(df50Lfdr) %>%
        mutate(df50Avg = cummean(df50Lfdr)) %>%
        # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
    powerRes$fdp50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 0)) /  length(which(totOut$df50Avg < 0.1))
    powerRes$power50df[sim_it] <- length(which(totOut$df50Avg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
    powerRes$nRej50df[sim_it] <- length(which(totOut$df50Avg < 0.1))
    } else {
      totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)
    }
  } else {totOut <- totOut %>% mutate(df50Lfdr = NA, df50Avg=NA)}
  
  # new method
  if (doNew) {
    initPiList <- list(c(0.82))
    for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
    initPiList[[8]] <- c(0.1)
    # the csmGmm package will add the appropriate 0s to initMuList
    initMuList <- list(matrix(data=0, nrow=3, ncol=1))
    for (i in 2:7) {
      initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
    }
    initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
    #newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))
    newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))
    
    # record
    totOut <- totOut %>% mutate(newLfdr = newRes$lfdrResults) %>%
      arrange(newLfdr) %>%
      mutate(newAvg = cummean(newLfdr)) %>%
      # don't forget to set it back to original index for further additions!!
      arrange(origIdx)
    powerRes$fdpNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$newAvg < 0.1))
    powerRes$powerNew[sim_it] <- length(which(totOut$newAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
    powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = newRes$lfdrResults))
    powerRes$nRejNew[sim_it] <- length(which(totOut$newAvg < 0.1))
  }
  cat('\n Done with ', sim_it, '\n')
  cat('\n \n \n \n \n \n \n \n \n \n \n')
}

setwd(outputDir)
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


