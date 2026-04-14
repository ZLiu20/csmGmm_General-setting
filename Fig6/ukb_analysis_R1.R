# Analysis of real data - pleiotropy and replication studies

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig6") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig6/ukb_analysis_R1.R")

# load libraries
library(dplyr)
library(data.table)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(usethis)
library(locfdr)
library(qch)
#library(adaFilter)
library(adaFilter, lib.loc = "/rsrch8/home/biostatistics/zliu20/R/x86_64-pc-linux-gnu-library")


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
define_H_space_R1 <- function(K,t,sameDirAlt=FALSE) {
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
  ## Revised on 2026/01/20
  if (sameDirAlt) {
    # H_annot[which(apply(abs(H_space), 1, sum) == K)] <- 1
    H_annot[which(
      apply(abs(H_space), 1, sum) > t &
        apply(abs(H_space), 1, sum) <= K & 
        apply(abs(H_space), 1, sum) == abs(apply(sign(H_space), 1, sum)) )] <- 1  
  } else { 
    H_annot[which(
      apply(abs(H_space), 1, sum) >t & apply(abs(H_space), 1, sum) <= K )] <- 1
  }
  
  return( list(H_space=H_space, H_annot=H_annot) )
}



################################################################################
emp_bayes_framework_R1 <- function(t_value, summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE, ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
  # Step 1 
  H_space_output <- define_H_space_R1(K = ncol(summary_tab), t = t_value, sameDirAlt = sameDirAlt)
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
#' Check for incongruous results in multiple dimensions
#'
#' @param zMatrix Matrix of test statistics
#' @param lfdrVec Vector of local false discovery rates
#'
#' @return Vector of indices of incongruous results
#'
#' @export
check_incongruous <- function(zMatrix, lfdrVec) {
  # Remove lfdr = 1
  lessThanOne <- which(lfdrVec < 0.99)
  if (length(lessThanOne) <= 1) {return(c())}
  
  zMatrix <- zMatrix[lessThanOne, ]
  lfdrVec <- lfdrVec[lessThanOne]
  
  # Do it in K^2 quadrants
  K <- ncol(zMatrix)
  quadrants <- expand.grid(rep(list(c(-1, 1)), K))
  
  badIdx <- c()
  for (quad_it in 1:nrow(quadrants)) {
    # Separate into quadrants
    idxVec <- 1:nrow(zMatrix)
    tempStats <- zMatrix
    tempLfdr <- lfdrVec
    
    for (k_it in 1:K) {
      if (class(tempStats)[1] == "numeric") {break}
      
      if (quadrants[quad_it, k_it] == -1) {
        toKeep <- which(tempStats[, k_it] < 0)
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      } else {
        toKeep <- which(tempStats[, k_it] > 0)
        idxVec <- idxVec[toKeep]
        tempLfdr <- tempLfdr[toKeep]
      }
      
      tempStats <- tempStats[toKeep, ]
    } # end finding quadrant
    
    if (length(idxVec) <= 1) {next}
    
    # Take absolute value
    tempStats <- abs(tempStats)
    
    # Prepare data
    tempDat <- tempStats %>% 
      as.data.frame(.data) %>%
      dplyr::mutate(lfdr = tempLfdr) %>%
      dplyr::mutate(idx = idxVec)
    
    # Check for incongruous results based on dimensions
    if (K == 2) {
      colnames(tempDat)[1:2] <- c("Z1", "Z2")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2))
      incongruousVec <- sapply(1:nrow(tempDat), 
                               FUN = find_2d, 
                               allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2)))
    } else if (K == 3) {
      colnames(tempDat)[1:3] <- c("Z1", "Z2", "Z3")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, dplyr::desc(.data$Z1), dplyr::desc(.data$Z2), dplyr::desc(.data$Z3))
      incongruousVec <- sapply(1:nrow(tempDat), 
                               FUN = find_3d, 
                               allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2, .data$Z3)))
    } else if (K == 4) {
      colnames(tempDat)[1:4] <- c("Z1", "Z2", "Z3", "Z4")
      tempDat <- tempDat %>%
        dplyr::arrange(tempLfdr, 
                       dplyr::desc(.data$Z1), 
                       dplyr::desc(.data$Z2), 
                       dplyr::desc(.data$Z3), 
                       dplyr::desc(.data$Z4))
      incongruousVec <- sapply(1:nrow(tempDat), 
                               FUN = find_4d, 
                               allTestStats = as.matrix(tempDat %>% dplyr::select(.data$Z1, .data$Z2, .data$Z3, .data$Z4)))
    } else {
      stop("only support for 2-4 dimensions right now")
    }
    
    # Get the bad indices
    badIdx <- c(badIdx, tempDat$idx[which(incongruousVec > 0)])
  }
  
  return(badIdx)
}


#'
find_2d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2]))
}

#'
find_3d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2] &
                 allTestStats[1:x, 3] < allTestStats[x, 3]))
}

#'
find_4d <- function(x, allTestStats) {
  length(which(
    allTestStats[1:x, 1] < allTestStats[x, 1] & 
      allTestStats[1:x, 2] < allTestStats[x, 2] & 
      allTestStats[1:x, 3] < allTestStats[x, 3] &
      allTestStats[1:x, 4] < allTestStats[x, 4]
  ))
}


################################################################################
run_qch_multi_alpha <- function(pmat, zmat, Hconfig_qch, H1config_qch, qch_test_name,
                                alpha_vec = c(0.01, 0.10)) {
  pmat <- as.matrix(pmat)
  zmat <- as.matrix(zmat)
  
  pmat[!is.finite(pmat)] <- 1
  pmat <- pmin(pmax(pmat, 0), 1)
  
  qch_fit <- qch::qch.fit(
    pValMat = as.data.frame(pmat),
    Hconfig = Hconfig_qch,
    copula = "gaussian",
    plotting = FALSE
  )
  
  out <- list()
  for (a in alpha_vec) {
    tag <- gsub("\\.", "", sprintf("%.2f", a))  # 0.01 -> "001", 0.10 -> "010"
    
    qch_res <- qch::qch.test(
      res.qch.fit = qch_fit,
      Hconfig = Hconfig_qch,
      Hconfig.H1 = H1config_qch,
      Alpha = a
    )
    
    lf <- qch_res$lFDR[[qch_test_name]]
    if (is.null(lf) || length(lf) != nrow(zmat)) {
      lf <- rep(NA_real_, nrow(zmat))
    } else {
      lf <- as.numeric(lf)
      lf[!is.finite(lf)] <- NA_real_
      lf <- pmin(pmax(lf, 0), 1)
    }
    
    incon <- if (all(is.na(lf))) NA_integer_ else {
      length(check_incongruous(zMatrix = zmat, lfdrVec = lf))
    }
    
    out[[paste0("lfdr_", tag)]] <- lf
    # out[[paste0("incon_", tag)]] <- incon
  }
  
  out
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
pleiotropy <- FALSE

doKernel <- FALSE
do50df <- FALSE
do7df <- FALSE
New <- FALSE
# doKernel <- TRUE
# do50df <- TRUE
# do7df <- TRUE
doQCH <- TRUE
doAdaFilter <- TRUE

t <- 1
if (aID == 1) {
  # three way - overall, cad, bmi
  pleiotropy <- TRUE
  replication <- FALSE
  cleanUKB <- fread(here::here(summaryStatDir, "bmi_with_overall.txt"))
  testDat <- cleanUKB %>% select(Zoverall, Zcad, Zbmi, pOverall, p_CAD, pBMI)
} else if (aID == 2) {
  # three way - overall, cad, ilcco
  replication <- TRUE
  pleiotropy <- FALSE
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


nZ <- ncol(testDat) / 2
zMat <- as.matrix(testDat[, 1:nZ, drop = FALSE])
pMat <- as.matrix(testDat[, (nZ + 1):ncol(testDat), drop = FALSE])

if (pleiotropy) {
  
  # kernel
  if (doKernel) {
    oldResKernel <- emp_bayes_framework_R1(
      t_value = t, summary_tab = zMat, sameDirAlt = replication,
      kernel = TRUE, joint = FALSE, ind = TRUE, dfFit = 7,
      Hdist_epsilon = 1e-2, checkpoint = TRUE
    )
    if (is.list(oldResKernel)) {
      kernelLfdr <- oldResKernel$lfdrVec
      inconKernel <- length(check_incongruous(zMatrix = zMat, lfdrVec = kernelLfdr))
    } else {
      kernelLfdr <- rep(NA_real_, nrow(zMat))
      inconKernel <- NA_integer_
    }
  }
  
  # 7df
  if (do7df) {
    oldRes7df <- emp_bayes_framework_R1(
      t_value = t, summary_tab = zMat, sameDirAlt = replication,
      kernel = FALSE, joint = FALSE, ind = TRUE, dfFit = 7,
      Hdist_epsilon = 1e-2, checkpoint = TRUE
    )
    if (is.list(oldRes7df)) {
      df7Lfdr <- oldRes7df$lfdrVec
      incon7df <- length(check_incongruous(zMatrix = zMat, lfdrVec = df7Lfdr))
    } else {
      df7Lfdr <- rep(NA_real_, nrow(zMat))
      incon7df <- NA_integer_
    }
  }
  
  # 50df
  if (do50df) {
    oldRes50df <- emp_bayes_framework_R1(
      t_value = t, summary_tab = zMat, sameDirAlt = replication,
      kernel = FALSE, joint = FALSE, ind = TRUE, dfFit = 50,
      Hdist_epsilon = 1e-2, checkpoint = TRUE
    )
    if (is.list(oldRes50df)) {
      df50Lfdr <- oldRes50df$lfdrVec
      incon50df <- length(check_incongruous(zMatrix = zMat, lfdrVec = df50Lfdr))
    } else {
      df50Lfdr <- rep(NA_real_, nrow(zMat))
      incon50df <- NA_integer_
    }
  }
  
  # New
  if (New) {
  initPiList <- list(c(0.82))
  for (i in 2:7) initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)
  initPiList[[8]] <- c(0.1)
  
  initMuList <- list(matrix(data = 0, nrow = 3, ncol = 1))
  for (i in 2:7) initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  initMuList[[8]] <- matrix(data = c(8, 8, 8), nrow = 3)
  
  newRes <- symm_fit_ind_EM_R1(
    t_value = t, testStats = zMat, initMuList = initMuList,
    initPiList = initPiList, sameDirAlt = replication, eps = newEps
  )
  inconnew <- if (is.list(newRes)) {
    length(check_incongruous(zMatrix = zMat, lfdrVec = newRes$lfdrResults))
  } else NA_integer_
  
  cat("[New method] DONE  ", format(Sys.time(), "%H:%M:%S"), "\n", sep = "")
  
  }
  
  
  # qch (two q values)
  
  if (doQCH) {
    Hconfig_qch <- qch::GetHconfig(ncol(zMat))
    H1config_qch <- qch::GetH1AtLeast(Hconfig_qch, AtLeast = 2)  # change if needed
    qch_test_name <- "AtLeast_2"
    
    qch_out <- tryCatch(
      run_qch_multi_alpha(
        pmat = pMat,
        zmat = zMat,
        Hconfig_qch = Hconfig_qch,
        H1config_qch = H1config_qch,
        qch_test_name = qch_test_name,
        alpha_vec = c(0.01, 0.10)
      ),
      error = function(e) list(
        lfdr_001 = rep(NA_real_, nrow(zMat)),
        incon_001 = NA_integer_,
        lfdr_010 = rep(NA_real_, nrow(zMat)),
        incon_010 = NA_integer_
      )
    )
    
    # save
    write.table(qch_out$lfdr_001, paste0(fnameRoot, "_qch01.txt"),
                append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(qch_out$incon_001, paste0(fnameRoot, "_inconqch01.txt"),
                append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)

    write.table(qch_out$lfdr_010, paste0(fnameRoot, "_qch10.txt"),
                append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(qch_out$incon_010, paste0(fnameRoot, "_inconqch10.txt"),
                append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
    
  
  # adaFilter
  if (doAdaFilter) {
    pmat_a <- pMat
    pmat_a[!is.finite(pmat_a)] <- 1
    pmat_a <- pmin(pmax(pmat_a, 0), 1)
    
    adaLfdr <- tryCatch({
      ao <- adaFilter::adaFilter(pmat_a, r = 2)
      if (!is.null(ao$adjusted.p) && length(ao$adjusted.p) == nrow(zMat)) as.numeric(ao$adjusted.p) else rep(NA_real_, nrow(zMat))
    }, error = function(e) rep(NA_real_, nrow(zMat)))
    
    adaLfdr[!is.finite(adaLfdr)] <- NA_real_
    adaLfdr <- pmin(pmax(adaLfdr, 0), 1)
    inconAda <- if (all(is.na(adaLfdr))) NA_integer_ else length(check_incongruous(zMat, adaLfdr))
  }
  
  
  if (doAdaFilter) {
    write.table(adaLfdr, paste0(fnameRoot, "_adaFilter.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(inconAda, paste0(fnameRoot, "_inconada.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  
  write.table(kernelLfdr, paste0(fnameRoot, "_kernel.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(inconKernel, paste0(fnameRoot, "_inconKernel.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(df7Lfdr, paste0(fnameRoot, "_df7.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(incon7df, paste0(fnameRoot, "_incon7df.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(df50Lfdr, paste0(fnameRoot, "_df50.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(incon50df, paste0(fnameRoot, "_incon50df.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(inconnew, paste0(fnameRoot, "_inconnew.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
}

if (replication) {
  
  if (doKernel) {
    oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = zMat, sameDirAlt = replication, kernel = TRUE, joint = FALSE, ind = TRUE, dfFit = 7, Hdist_epsilon = 1e-2, checkpoint = TRUE)
    if (is.list(oldResKernel)) {
      kernelLfdr <- oldResKernel$lfdrVec
      inconKernel <- length(check_incongruous(zMat, kernelLfdr))
    } else {
      kernelLfdr <- rep(NA_real_, nrow(zMat)); inconKernel <- NA_integer_
    }
  }
  
  if (do7df) {
    oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = zMat, sameDirAlt = replication, kernel = FALSE, joint = FALSE, ind = TRUE, dfFit = 7, Hdist_epsilon = 1e-2, checkpoint = TRUE)
    if (is.list(oldRes7df)) {
      df7Lfdr <- oldRes7df$lfdrVec
      incon7df <- length(check_incongruous(zMat, df7Lfdr))
    } else {
      df7Lfdr <- rep(NA_real_, nrow(zMat)); incon7df <- NA_integer_
    }
  }
  
  if (do50df) {
    oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = zMat, sameDirAlt = replication, kernel = FALSE, joint = FALSE, ind = TRUE, dfFit = 50, Hdist_epsilon = 1e-2, checkpoint = TRUE)
    if (is.list(oldRes50df)) {
      df50Lfdr <- oldRes50df$lfdrVec
      incon50df <- length(check_incongruous(zMat, df50Lfdr))
    } else {
      df50Lfdr <- rep(NA_real_, nrow(zMat)); incon50df <- NA_integer_
    }
  }
  
  initPiList <- list(c(0.82))
  for (i in 2:7) initPiList[[i]] <- c(0.08 / 12, 0.08 / 12)
  initPiList[[8]] <- c(0.1)
  initMuList <- list(matrix(data = 0, nrow = 3, ncol = 1))
  for (i in 2:7) initMuList[[i]] <- cbind(rep(2, 3), rep(5, 3))
  initMuList[[8]] <- matrix(data = c(8, 8, 8), nrow = 3)
  
  newRes <- symm_fit_ind_EM_R1(
    t_value = t, testStats = zMat, initMuList = initMuList,
    initPiList = initPiList, sameDirAlt = replication, eps = newEps
  )
  inconnew <- if (is.list(newRes)) length(check_incongruous(zMat, newRes$lfdrResults)) else NA_integer_
  
  write.table(kernelLfdr, paste0(fnameRoot, "_kernel.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(inconKernel, paste0(fnameRoot, "_inconKernel.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(df7Lfdr, paste0(fnameRoot, "_df7.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(incon7df, paste0(fnameRoot, "_incon7df.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(df50Lfdr, paste0(fnameRoot, "_df50.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(incon50df, paste0(fnameRoot, "_incon50df.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(newRes$lfdrResults, paste0(fnameRoot, "_newlfdr.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(inconnew, paste0(fnameRoot, "_inconnew.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(do.call(cbind, newRes$muInfo), paste0(fnameRoot, "_muInfo.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(do.call(cbind, newRes$piInfo), paste0(fnameRoot, "_piInfo.txt"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
}





