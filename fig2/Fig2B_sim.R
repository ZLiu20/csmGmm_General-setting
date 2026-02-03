# For Figure 2B

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1C_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig2/Fig2B_sim.R")

# # load libraries
# library(mvtnorm)
# library(data.table)
# library(bindata)
# library(tidyverse)
# library(devtools)
# library(ks)
# library(csmGmm)
# library(here)
# library(locfdr)


# load libraries
library(dplyr)
library(mvtnorm)
library(data.table)
library(bindata)
library(magrittr)
library(rje)
library(ks)
library(csmGmm)
library(here)
library(locfdr)



# define the  function

################################################################################
find_max_means_R1 <- function(muInfo) {
  
  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  S1 <- c(1,2,4,5,5,8,9,10)+1
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
  H_annot[which(
    t+1<=apply(abs(H_space), 1, sum) & 
      apply(abs(H_space), 1, sum) <= K &
      apply(H_space, 1, function(row) {
        (row[1]*row[2]) +(row[3]*row[4])
      })!=0
  )] <- 1
  
  
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
                            apply(H_space, 1, function(row) {
                              (row[1]*row[2]) +(row[3]*row[4])
                            })==0 )) 
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
  sum_12 <- apply(Hmat[, 1:K], 1, function(row) row[1] * row[2])
  sum_34 <- apply(Hmat[, 1:K], 1, function(row) row[3] * row[4])
  
  symAltVec <- ifelse((t+1 <= apply(abs(Hmat[, 1:K]), 1, sum)) & 
                        (apply(abs(Hmat[, 1:K]), 1, sum) <= K) &  
                        (sum_12 + sum_34 !=0) , 1, 0)
  
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
        S2 <- c(3,7,11,12,13,14,15)
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

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Fig2", "output")
outName <- paste0(outputDir, "/Fig2B_aID", aID, ".txt")

# option to save or load intermediate data to save time
loadData <- FALSE
saveData <- FALSE
# these names are for if saveData <- TRUE
testStatsName <- here::here(outputDir, "Fig2B_allZ")
betaName <- here::here(outputDir, "Fig2B_allBeta")

# simulation parameters start here
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
t <- 1
nSets <- nSNPs / setSize
nSims <- 10
margprob <- rep(0.3, setSize)
simsPerEffSize <- 40
effSizeMult <- ceiling(aID / simsPerEffSize)

betaMin <- c(0.14, 0.14, 0.14, 0.14)
betaMax <- c(0.18, 0.18, 0.18, 0.18)
                       
# betaMin <- c(0.14, 0.18, 0.14, 0.18)
# betaMax <- c(0.14, 0.18, 0.14, 0.18)
                       
beta0 <- -1

# determines how many signals there are
# pi1 <- 0.01 * effSizeMult
# pi11 <- 3 * pi1^2
# pi111 <- pi11 / 2
# pi00 <- 1 - 3 * pi1 - 3 * pi11 - pi111
# sProp <- c(pi00, 3 * pi1, 3 * pi11, pi111)

pi1 <- 0.01 * effSizeMult
pi11 <- 3 * pi1^2
pi111 <- pi11 / 2
pi1111 <- pi111 / 2  
pi00 <- 1 - 3 * pi1 - 3 * pi11 - pi111 - pi1111  
sProp <- c(pi00, 3 * pi1, 3 * pi11, pi111, pi1111) 

hMat <- expand.grid(c(-1, 0, 1), c(-1, 0, 1)) %>%
  as.data.frame(.) %>%
  mutate(s = abs(Var1) + abs(Var2)) %>%
  arrange(s)
number <- c()
for (s_it in 0:max(hMat$s)) {
  numRows <- length(which(hMat$s == s_it))
  number <- c(number, rep(sProp[s_it + 1] * nSNPs / numRows, numRows))
}
hMat <- hMat %>% mutate(number = number)

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims),  minEff1=pi1,
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA, pi0cTrue=NA,
                       nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA, nRejNew=NA,
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
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else {
    # hold test statistics and signals
    allZ <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    allBeta <- matrix(data=NA, nrow=nSNPs, ncol=nDims)
    
    # select signal locations
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims = nDims)
    
    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets)  {
      
      statsMat <- matrix(data=NA, nrow=setSize, ncol=nDims)
      
      # generate coefficient matrix
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize, 
                          betaMin=betaMin, betaMax=betaMax)
      
      # two dimensional mediation case data generation for k=1 dimension
      tempG <- sapply(X=margprob, FUN=rbinom, n=n, size=2)
      adjG <- sweep(tempG, MARGIN=2, STATS=apply(tempG, 2, mean), FUN="-")
      
      ### M1
      tempAlpha1 <- coefMat[, 1] 
      tempM1 <- sweep(tempG, MARGIN=2, STATS=tempAlpha1, FUN="*") + matrix(data=rnorm(n=n * nrow(coefMat)), nrow=n, ncol=nrow(coefMat))
      adjM1 <- sweep(tempM1, MARGIN = 2, STATS = apply(tempM1, 2, mean), FUN = "-")
      sigSqHat1 <- apply(adjM1, 2, myvar_fun)
      
      # calculate test statistics for k=1
      tempNum1 <- apply(adjG * adjM1, 2, sum)
      tempDenom1 <- sqrt(apply(adjG^2, 2, sum) * sigSqHat1)
      statsMat[, 1] <- tempNum1 / tempDenom1
      
      # mediation data generation for k=2 dimension 
      tempBeta1 <- coefMat[, 2]
      tempEta1 <- sweep(tempM1, MARGIN=2, STATS=tempBeta1, FUN="*") + matrix(data=beta0, nrow=nrow(tempM1), ncol=ncol(tempM1))
      tempMu1 <- rje::expit(as.numeric(tempEta1))
      tempY1 <- rbinom(n=length(tempMu1), size=1, prob=tempMu1)
      yMat1 <- matrix(data=tempY1, nrow=nrow(tempEta1), ncol=ncol(tempEta1), byrow=FALSE)
      
      # calculate test statistics for k=2 
      for (test_it in 1:ncol(yMat1)) {
        tempMod1 <- glm(yMat1[, test_it] ~ tempG[, test_it] + tempM1[, test_it], family=binomial)
        statsMat[test_it, 2] <- summary(tempMod1)$coefficients[3, 3]
      } 
      
      ### M2
      tempAlpha2 <- coefMat[, 3]
      tempM2 <- sweep(tempG, MARGIN=2, STATS=tempAlpha2, FUN="*") + matrix(data=rnorm(n=n * nrow(coefMat)), nrow=n, ncol=nrow(coefMat))
      adjM2 <- sweep(tempM2, MARGIN = 2, STATS = apply(tempM2, 2, mean), FUN = "-")
      sigSqHat2 <- apply(adjM2, 2, myvar_fun)
      
      # calculate test statistics for k=1
      tempNum2 <- apply(adjG * adjM2, 2, sum)
      tempDenom2 <- sqrt(apply(adjG^2, 2, sum) * sigSqHat2)
      statsMat[, 3] <- tempNum2 / tempDenom2
      
      # mediation data generation for k=2 dimension 
      tempBeta2 <- coefMat[, 4]
      tempEta2 <- sweep(tempM2, MARGIN=2, STATS=tempBeta2, FUN="*") + matrix(data=beta0, nrow=nrow(tempM2), ncol=ncol(tempM2))
      tempMu2 <- rje::expit(as.numeric(tempEta2))
      tempY2 <- rbinom(n=length(tempMu2), size=1, prob=tempMu2)
      yMat2 <- matrix(data=tempY2, nrow=nrow(tempEta2), ncol=ncol(tempEta2), byrow=FALSE)
      
      
      # calculate test statistics for k=2 
      for (test_it in 1:ncol(yMat2)) {
        tempMod2 <- glm(yMat2[, test_it] ~ tempG[, test_it] + tempM2[, test_it], family=binomial)
        statsMat[test_it, 4] <- summary(tempMod2)$coefficients[3, 3]
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
      write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    } 
  }
  
  # number of signals and causal SNPs
  #causalVec <- as.numeric(allBeta[, 1] != 0 & allBeta[, 2] != 0)
  causalVec <- as.numeric(
    rowSums(allBeta != 0) >= t+1 & rowSums(allBeta != 0) <= nDims &
      ( allBeta[, 1] * allBeta[, 2] + allBeta[, 3] * allBeta[, 4] != 0 ) )
  
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
    #                                     dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, sameDirAlt=TRUE, kernel = TRUE, joint=FALSE, ind = TRUE,
                                           dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 0)) / length(which(totOut$kernelAvg < qvalue))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue & totOut$causal == 1)) / sum(causalVec)
      # powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$inconKernel[sim_it] <- length(check_incongruous_R1(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec,t_value=t))
      powerRes$nRejKernel[sim_it] <- length(which(totOut$kernelAvg < qvalue)) 
    } else {
      totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)
    }
  } else {totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)}
  
  # old method - 7 df
  if (do7df) { 
    # oldRes7df <- emp_bayes_framework(summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE, 
    #                                  dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, sameDirAlt=TRUE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                        dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
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
    #                                   dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
    oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, sameDirAlt=TRUE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                         dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
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
    # initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
    # initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
    #                    matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))
    
    initPiList <- list(c(0.82))
    for (i in 2:(2^nDims-1)) {initPiList[[i]] <- c(0.08 / 28, 0.08 / 28)}
    initPiList[[2^nDims]] <- c(0.1)
    # the symm_fit_ind.R code will add the appropriate 0s to initMuList
    initMuList <- list(matrix(data=0, nrow=nDims, ncol=1))
    for (i in 2:(2^nDims-1)) {
      initMuList[[i]] <- cbind(rep(2, nDims), rep(5, nDims))
    }
    initMuList[[2^nDims]] <- matrix(data=c(8, 8, 8, 8), nrow=nDims)
    
    # initMuList <- list(matrix(data=rep(0, nDims), nrow=nDims, ncol=1))
    # for (i in 2:(2^nDims)) {
    #   initMuList[[i]] <- matrix(data=rep(3, nDims), nrow=nDims, ncol=1)
    # }
    
    # newRes <- symm_fit_ind_EM(testStats = allZ, initMuList = initMuList, initPiList = initPiList, eps=10^(-5))
    newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = allZ, initMuList = initMuList,
                                 initPiList = initPiList, sameDirAlt=TRUE, eps=10^(-3))
    
    # record
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
  }
  cat('\n Done with ', sim_it, '\n')
}

# save
write.table(powerRes, outName, append=F, quote=F, row.names=F, col.names=T, sep='\t')


