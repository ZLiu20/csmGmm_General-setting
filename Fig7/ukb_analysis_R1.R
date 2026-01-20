# Analysis of real data - pleiotropy and replication studies

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig1/Fig1A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig4/ukb_analysis_R1.R")

# load libraries
library(dplyr)
library(data.table)
library(devtools)
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
  S1 <- c(1,2,4)+1
  # just keep finding the max
  maxMeans <- rep(0, K) 
  for (element_it in S1) {
    tempMat <- cbind(muInfo[[element_it]], maxMeans)
    maxMeans <- apply(tempMat, 1, max)
  }
  # return K*1 vector
  return(maxMeans)
}


################################################################################
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
    apply(abs(H_space), 1, sum) > t &
      apply(abs(H_space), 1, sum) <= K & 
      apply(H_space, 1, sum) == abs(apply(sign(H_space), 1, sum))
  )] <- 1  
  
  return( list(H_space=H_space, H_annot=H_annot) )
}



################################################################################
emp_bayes_framework_R1 <- function(t_value, summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE, ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
  # Step 1 
  H_space_output <- define_H_space_R1(K = ncol(summary_tab), t = t_value)
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
  locfdr_output <- locfdr_estim(summary_tab=summary_tab, kernel=kernel, joint=joint, ind=ind, dfFit = dfFit, nulltype=nulltype)
  
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
symm_fit_ind_EM_R1 <- function(t_value, testStats, initMuList, initPiList, sameDirAlt=TRUE, eps = 10^(-5), checkpoint=TRUE) {
  
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
    slVec <- slVec + abs(Hmat[, k_it])
  }
  # symmetric alternative
  sum_Hmat <- apply(abs(Hmat[, 1:K]), 1, sum)
  symAltVec <- ifelse(
    sum_Hmat > t & 
      sum_Hmat <= K & 
      sum_Hmat == abs(apply(sign(Hmat), 1, sum)), 1, 0)  # revised on 15/11/2025
  
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
    # nullCols <- which(allPi[, 3] < K)
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
outputDir <- here::here("Fig4", "output")
fnameRoot <- paste0(outputDir, "/Fig4_data_aID", aID)

# where is the data
summaryStatDir <- here::here("Data")

# controls convergence of EM algorithms
oldEps <- 0.01
newEps <- 10^(-3)

# controls which analysis to perform
replication <- FALSE
threeway <- FALSE

doKernel <- TRUE
do50df <- TRUE
do7df <- TRUE

t <- 1
if (aID == 1) {
  # three way - overall, cad, bmi
  threeway <- TRUE
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zcad, Zbmi, pOverall, p_CAD, pBMI)
} else if (aID == 2) {
  # three way - overall, cad, ilcco
  replication <- TRUE
  cleanUKB <- fread(here::here(summaryStatDir, "lc_overall_three.txt"))
  testDat <- cleanUKB %>% select(zILCCO, zUKB, zMVP, pILCCO, pUKB, pMVP)
} 

testDat <- testDat %>% as.matrix(.)

# adjust so test statistics are not too large for R
testDat <- testDat %>% as.matrix(.)
for (col_it in 1:(ncol(testDat)/2)) {
  tooBig <- which(testDat[, col_it] > 8.1)
  tooSmall <- which(testDat[, col_it] < -8.1)
  if (length(tooBig) > 0) {
    testDat[tooBig, col_it] <- 8.1
  }
  if (length(tooSmall) > 0) {
    testDat[tooSmall, col_it] <- -8.1
  }
}
# adjust so p-values are not 0 or 1
for (col_it in ((ncol(testDat)/2)+1):ncol(testDat)) {
  tooBig <- which(testDat[, col_it] == 1)
  tooSmall <- which(testDat[, col_it] == 0)
  minVal <- min(testDat[which(testDat[, col_it] > 0), col_it]) 
  if (length(tooBig) > 0) {                
    testDat[tooBig, col_it] <- 0.999        
  }        
  if (length(tooSmall) > 0) {                
    testDat[tooSmall, col_it] <- minVal          
  }
}


# run kernel
if (doKernel) {
  oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = TRUE, joint=FALSE, ind = TRUE,
                                         dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  # oldResKernel <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = TRUE, joint=FALSE, ind = TRUE,
  #                                   dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldResKernel)[1] != "list") {
    kernelLfdr <- rep(NA, nrow(testDat))
  } else {
    kernelLfdr <- oldResKernel$lfdrVec
  }
  # save
  write.table(kernelLfdr, paste0(fnameRoot, "_kernel.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# run 7 df
if (do7df) {
  # oldRes7df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
  #                                  dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
                                      dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes7df)[1] != "list") {
    df7Lfdr <- rep(NA, nrow(testDat))
  } else {
    df7Lfdr <- oldRes7df$lfdrVec
  }
  # save
  write.table(df7Lfdr, paste0(fnameRoot, "_df7.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# run 50 df
if (do50df) {
  # oldRes50df <- emp_bayes_framework(summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
  #                                 dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = testDat[, 1:(ncol(testDat)/2)], sameDirAlt = replication, kernel = FALSE, joint=FALSE, ind = TRUE,
                                       dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  if (class(oldRes50df)[1] != "list") {
    df50Lfdr <- rep(NA, nrow(testDat))
  } else {
    df50Lfdr <- oldRes50df$lfdrVec
  }
  # save
  write.table(df50Lfdr, paste0(fnameRoot, "_df50.txt"), append=F, quote=F, row.names=F, col.names=T)
}

# 3D cases pleiotropy
if (threeway) {
  initPiList <- list(c(0.82))
  for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiList[[8]] <- c(0.1)
  tempH <- expand.grid(c(0, 1), c(0, 1), c(0, 1)) %>%
    mutate(s = Var1 + Var2 + Var3) %>%
    arrange(s) %>%
    select(-s) %>%
    as.matrix(.)
  initMuList <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  }
  initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
  
  newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = testDat[, 1:3], initMuList = initMuList, 
                               initPiList = initPiList, eps = newEps)
}


# 3D cases replication
if (replication) {
  initPiList <- list(c(0.82))
  for (i in 2:7) {initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiList[[8]] <- c(0.1)
  # the csmGmm package will add the appropriate 0s to initMuList
  initMuList <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  }
  initMuList[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
  
  newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = testDat[, 1:3], initMuList = initMuList, 
                               initPiList = initPiList, sameDirAlt=TRUE, eps = newEps)
  
}

# save
write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)







