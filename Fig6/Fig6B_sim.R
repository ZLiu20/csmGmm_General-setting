# Figure 6B

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig6/Fig6A_sim.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig6/Fig6B_sim.R")

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
library(Matrix)    
#library(corpcor)   


## Define the supporting functions (unchanged)

################################################################################
find_max_means_R1 <- function(muInfo) {
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  S1 <- c(1,2,4)+1
  maxMeans <- rep(0, K) 
  for (element_it in S1) {
    tempMat <- cbind(muInfo[[element_it]], maxMeans)
    maxMeans <- apply(tempMat, 1, max)
  }
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
  H_annot[which(t+1 <= apply(abs(H_space), 1, sum) & apply(abs(H_space), 1, sum) <= K)] <- 1
  return( list(H_space=H_space, H_annot=H_annot) )
}

################################################################################
emp_bayes_framework_R1 <- function(t_value, summary_tab, sameDirAlt=FALSE, nulltype=1, kernel=TRUE, joint=FALSE,
                                   ind=TRUE, dfFit=7, Hdist_epsilon=10^(-2), checkpoint=TRUE) {
  H_space_output <- define_H_space_R1(K = ncol(summary_tab), t=t_value)
  H_space <- as.matrix(H_space_output$H_space)
  H_annot <- H_space_output$H_annot
  if (sameDirAlt) {
    null_rows <- which(apply(abs(H_space), 1, sum) <= t_value |
                         (apply(abs(H_space), 1, sum) > t_value &
                            apply(abs(H_space), 1, sum) <= ncol(summary_tab) &
                            apply(abs(H_space), 1, sum) != abs(apply(sign(H_space), 1, sum))))
  } else { 
    null_rows <- which(H_annot == 0)
  }
  
  locfdr_output <- locfdr_estim(summary_tab=summary_tab, kernel=kernel, joint=joint, ind=ind, dfFit = dfFit, nulltype = nulltype)
  pi0_estim <- locfdr_output$pi0_estim
  if (length(which(pi0_estim >= 1 | pi0_estim <= 0)) > 0) {
    pi0_estim[which(pi0_estim >= 1 | pi0_estim <= 0)] <- 0.99
  }
  marg_pmf_tab <- locfdr_output$marg_pmf_tab
  if (is.null(marg_pmf_tab)) {return(-1)} 
  if (ncol(marg_pmf_tab) < 2*ncol(summary_tab)) {
    return(-1)
  }
  
  H_dist <- calc_H_dist_indep(pi0_estim=pi0_estim, H_space=H_space)
  cond_pmfs <- calc_cond_pmfs(summary_tab=summary_tab, marg_pmf_tab=marg_pmf_tab, pi0_estim=pi0_estim)
  null_pmf_tab <- cond_pmfs$null_pmf_tab
  neg_pmf_tab <- cond_pmfs$neg_pmf_tab
  pos_pmf_tab <- cond_pmfs$pos_pmf_tab
  
  if (is.na(sum(neg_pmf_tab)) | is.na(sum(neg_pmf_tab))) {return(-1)}
  
  binned_tab <- cond_pmfs$binned_tab
  EMoutput <- run_EM_pmf(binned_dat=binned_tab, H_dist=H_dist, H_space=H_space,
                         null_pmf_tab=null_pmf_tab, neg_pmf_tab=neg_pmf_tab,
                         pos_pmf_tab=pos_pmf_tab, pi0_estim=pi0_estim, epsilon=Hdist_epsilon)
  
  loc_fdr_num <- apply(EMoutput$jointMat[, null_rows], 1, sum)
  lfdrVec <- loc_fdr_num / EMoutput$probZ
  return(list(marg_pmf_tab = marg_pmf_tab, null_pmf_tab = null_pmf_tab, pi0_estim = pi0_estim, Hdist_final = EMoutput$H_dist,
              neg_pmf_tab = neg_pmf_tab, pos_pmf_tab = pos_pmf_tab, lfdrVec = lfdrVec))
}

################################################################################
symm_fit_cor_EM_R1 <- function(t_value, testStats, corMat, initMuList, initPiList, eps = 10^(-5), checkpoint=TRUE) {
  sigInv = solve(corMat)
  J <- nrow(testStats)
  K <- ncol(testStats)
  B <- 2^K - 1
  L <- 3^K - 1
  t <- t_value
  Hmat <- expand.grid(rep(list(-1:1), K))
  blVec <- slVec <- rep(0, nrow(Hmat))
  for (k_it in K:1) {
    blVec <- blVec + 2^(K - k_it) * abs(Hmat[, k_it])
    slVec <- slVec + abs(Hmat[, k_it])
  }
  Hmat <- Hmat %>% dplyr::mutate(bl = blVec) %>%
    dplyr::mutate(sl = slVec) %>%
    dplyr::arrange(.data$bl, .data$Var1, .data$Var2) %>%
    dplyr::mutate(l = 0:(nrow(.) - 1))
  
  # initialize
  # the muInfo and piInfo are lists that hold the information in a more compact manner.
  # the allMu and allPi are matrices that repeat the information so the calculations can be performed faster.
  muInfo <- initMuList
  piInfo <- initPiList
  oldParams <- c(unlist(piInfo), unlist(muInfo))
  MbVec <- sapply(piInfo, FUN=length)
  
  diffParams <- 10
  iter <- 0
  while (diffParams > eps) {
    allPi <- c(piInfo[[1]], 0, 0, 1, 0)
    allMu <- rep(0, K)
    
    for (b_it in 1:B) {
      tempH <- Hmat %>% filter(.data$bl == b_it)
      for (m_it in 1:MbVec[b_it + 1]) {
        allPi <- rbind(allPi, cbind(rep(piInfo[[b_it + 1]][m_it] / (2^tempH$sl[1]), nrow(tempH)), tempH$bl,
                                    tempH$sl, rep(m_it, nrow(tempH)), tempH$l))
        for (h_it in 1:nrow(tempH)) {
          allMu <- cbind(allMu, unlist(tempH %>% dplyr::select(-.data$bl, -.data$sl, -.data$l) %>%
                                         dplyr::slice(h_it)) * muInfo[[b_it + 1]][, m_it])
        }
      }
    }
    colnames(allPi) <- c("Prob", "bl", "sl", "m", "l")
    
    conditionalMat <- sapply(X=data.frame(allMu), FUN=calc_dens_cor, Zmat = testStats, corMat = corMat) %>%
      sweep(., MARGIN=2, STATS=allPi[, 1], FUN="*")
    
    probZ <- apply(conditionalMat, 1, sum)
    AikMat <- conditionalMat %>% sweep(., MARGIN=1, STATS=probZ, FUN="/")
    Aik_alln <- apply(AikMat, 2, sum) / J
    
    for (b_it in 0:B) {
      for (m_it in 1:MbVec[b_it + 1]) {
        tempIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        piInfo[[b_it + 1]][m_it] <- sum(Aik_alln[tempIdx])
      }
    }
    
    for (b_it in 1:B) {
      tempHmat <- Hmat %>% filter(.data$bl == b_it)
      for (m_it in 1:MbVec[b_it + 1]) {
        tempRightSum <- rep(0, K)
        tempLeftSum <- rep(0, K)
        AikIdx <- which(allPi[, 2] == b_it & allPi[, 4] == m_it)
        for (idx_it in 1:length(AikIdx)) {
          tempAik <- AikIdx[idx_it]
          tempHvec <- tempHmat %>% dplyr::select(-.data$bl, -.data$sl, -.data$l) %>%
            dplyr::slice(idx_it) %>% unlist(.)
          LsigInvL <- diag(tempHvec) %*% sigInv %*% diag(tempHvec)
          
          tempLeftSum <- tempLeftSum + colSums(AikMat[, tempAik] * sweep(x = testStats %*% sigInv, MARGIN = 2,
                                                                         STATS = tempHvec, FUN="*"))
          tempRightSum <- tempRightSum + (J * Aik_alln[tempAik]) * diag(LsigInvL)
        }
        
        if (b_it == B) {
          muInfo[[b_it + 1]][, m_it] <- solve(diag(tempRightSum)) %*% tempLeftSum
        } else {
          whichZero <- which(tempRightSum == 0)
          tempRightSum[whichZero] <- 1
          muInfo[[b_it + 1]][, m_it] <- tempLeftSum / tempRightSum
        }
        
        S2 <- c(3,5,6,7)
        if (b_it %in% S2) {
          maxMeans <- find_max_means_R1(muInfo)
          whichSmaller <- which(muInfo[[b_it+1]][, m_it] < maxMeans)
          if (length(whichSmaller) > 0) {
            muInfo[[b_it+1]][whichSmaller, m_it] <- maxMeans[whichSmaller]
          }
        }
      }
    }
    
    allParams <- c(unlist(piInfo), unlist(muInfo))
    diffParams <- sum((allParams - oldParams)^2)
    oldParams <- allParams
    iter <- iter + 1
    if (checkpoint) {
      cat(iter, " - ", diffParams, "\n")
    }
  }
  
  nullCols <- which(allPi[, 3] <= t)
  probNull <- apply(conditionalMat[, nullCols], 1, sum)
  lfdrResults <- probNull / probZ
  
  return(list(piInfo = piInfo, muInfo = muInfo, iter = iter, lfdrResults = lfdrResults))
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

find_2d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2]))
}


find_3d <- function(x, allTestStats) {
  length(which(allTestStats[1:x, 1] < allTestStats[x, 1] & allTestStats[1:x, 2] < allTestStats[x, 2] &
                 allTestStats[1:x, 3] < allTestStats[x, 3]))
}


################################################################################
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

outputDir <- here::here("Fig6", "output")
outName <- paste0(outputDir, "/Fig6B_aID", aID, ".txt")

loadData <- FALSE
saveData <- FALSE
testStatsName <- here::here(outputDir, "Fig6B_allZ")
betaName <- here::here(outputDir, "Fig6B_allBeta")

# Simulation parameters
outcomeCor <- 0.1
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
nDims <- 3   # Now 3 dimensions
t <- 1
nSets <- nSNPs / setSize
nSims <- 5
margprob <- rep(0.3, setSize)
simsPerEffSize <- 20
effSizeMult <- ceiling(aID / simsPerEffSize)
betaMin <- rep(0.14 + 0.01 * effSizeMult, nDims)
betaMax <- rep(0.24 + 0.01 * effSizeMult, nDims)
beta0 <- -1

# Signal proportion for s = 0,1,2,3
sProp <- c(0.9582, 0.04, 0.0012, 0.0006)
# sProp <- c(0.9579, 0.04, 0.0012, 0.0006, 0.0003)


# ---- MODIFIED: 3D hMat ----
hMat <- expand.grid(rep(list(c(-1, 0, 1)), nDims)) %>%
  as.data.frame() %>%
  setNames(paste0("Var", 1:nDims)) %>%
  mutate(s = rowSums(abs(select(., starts_with("Var"))))) %>%
  arrange(s)

number <- c()
for (s_it in 0:max(hMat$s)) {
  numRows <- sum(hMat$s == s_it)
  number <- c(number, rep(sProp[s_it + 1] * nSNPs / numRows, numRows))
}
hMat$number <- number

# tells rmvbin how to generate correlated binary outcomes
sigmaValsMat <- matrix(outcomeCor, nrow = nDims, ncol = nDims)
diag(sigmaValsMat) <- 1


# generate two correlated outcomes for one subject
gen_one_subject <- function(x, sigmaValsMat) {
  # x should be rounded
  return(rmvbin(n=1, margprob=x, sigma=sigmaValsMat))
}

# record results here
powerRes <- data.frame(nCausal=rep(NA, nSims), minEff1=betaMin[1],
                       seed=NA, pi0aTrue=NA, pi0bTrue=NA, pi0cTrue = NA,
                       nRejDACT=NA, nRejHDMT=NA, nRejKernel=NA, nRej7df=NA, nRej50df=NA, nRejNew=NA,
                       powerDACT=NA, powerHDMT=NA, powerKernel=NA, power7df=NA,
                       power50df=NA, powerNew=NA, fdpDACT=NA, fdpHDMT=NA, fdpKernel=NA, fdp7df=NA, fdp50df=NA, fdpNew=NA,
                       inconKernel=NA, incon7df=NA, incon50df=NA, inconNew=NA)

for (sim_it in 1:nSims) {
  set.seed(aID * 10^5 + sim_it)
  powerRes$seed[sim_it] <- aID * 10^5 + sim_it
  
  if (loadData) {
    setwd(outputDir)
    allZ <- fread(paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
    allBeta <- fread(paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), data.table=F)
  } else {
    allZ <- matrix(NA, nrow = nSNPs, ncol = nDims)
    allBeta <- matrix(NA, nrow = nSNPs, ncol = nDims)
    
    sigLocsMat <- allocate_sigs(hMat = hMat, nSNPs = nSNPs, nDims = nDims)
    
    # generate data in multiple sets - faster than all at once
    for (set_it in 1:nSets) {
      statsMat <- matrix(NA, nrow = setSize, ncol = nDims)
      coefMat <- set_beta(sigLocsMat = sigLocsMat, set_it = set_it, setSize = setSize,
                          betaMin = betaMin, betaMax = betaMax)
      
      tempG <- sapply(X = margprob, FUN = rbinom, n = n, size = 2)
      etaMat <- tempG %*% coefMat + matrix(beta0, nrow = n, ncol = ncol(coefMat))
      muMat <- rje::expit(etaMat)
      
      # make the outcome
      roundedMu <- round(muMat, digits=2)
      for (dim_it in 1:nDims) {
        tooSmall <- which(roundedMu[, dim_it] == 0)
        if (length(tooSmall) > 0) {
          roundedMu[tooSmall, dim_it] <- 0.01
        }
        tooBig <- which(roundedMu[, dim_it] == 1)
        if (length(tooBig) > 0) {
          roundedMu[tooBig, dim_it] <- 0.99
        }
      }
      
      yMat <- t(apply(muMat, 1, gen_one_subject, sigmaValsMat = sigmaValsMat))
      
      # calculate test statistics
      for (dim_it in 1:nDims) {
        for (snp_it in 1:setSize) {
          tempMod <- glm(yMat[, dim_it] ~ tempG[, snp_it], family = binomial)
          statsMat[snp_it, dim_it] <- summary(tempMod)$coefficients[2, 3]
        }
      }
      
      startIdx <- (set_it - 1) * setSize + 1
      endIdx <- set_it * setSize
      allZ[startIdx:endIdx, ] <- statsMat
      allBeta[startIdx:endIdx, ] <- coefMat
      
      if(set_it %% 10 == 0) cat("Set", set_it, "\n\n\n")
    }
    
    if (saveData) {
      setwd(outputDir)
      write.table(allZ, paste0(testStatsName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
      write.table(allBeta, paste0(betaName, "_aID", aID, "_sim", sim_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
    }
  }
  
  # Causal definition: at least t+1 = 2 non-zero effects
  causalVec <- as.numeric(rowSums(allBeta != 0) >= (t + 1) & rowSums(allBeta != 0) <= nDims)
  powerRes$nCausal[sim_it] <- sum(causalVec)
  powerRes$pi0aTrue[sim_it] <- sum(allBeta[, 1] != 0)
  powerRes$pi0bTrue[sim_it] <- sum(allBeta[, 2] != 0)
  powerRes$pi0cTrue[sim_it] <- sum(allBeta[, 3] != 0)
  
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
      powerRes$fdpHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1 & totOut$causal == 0)) /  length(which(totOut$hdmtRes < 0.1))
      powerRes$powerHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$nRejHDMT[sim_it] <- length(which(totOut$hdmtRes < 0.1))
      
    } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  } else {totOut <- totOut %>% mutate(hdmtRes = NA)}
  
  # DACT
  if (class(nullprop)[1] != "list") {
    nullprop <- NULL
  }
  if (doDACT) {
    DACTout <- tryCatch(DACT_noEst(p_a = allP[, 1], p_b = allP[, 2], nullEst = nullprop, correction="JC"),
                        error = function(e) e, warning=function(w) w)
    if (class(DACTout)[1] %in% c("simpleWarning", "simpleError")) {
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
  
  
  # Kernel, 7df, 50df, New method — all support nDims
  if (doKernel) {
    oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = TRUE, joint=FALSE, ind = TRUE,
                                           dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=FALSE)
    if (class(oldResKernel)[1] == "list") {
      totOut <- totOut %>% mutate(kernelLfdr = oldResKernel$lfdrVec) %>%
        arrange(kernelLfdr) %>%
        mutate(kernelAvg = cummean(kernelLfdr)) %>% # don't forget to set it back to original index for further additions!!
        arrange(origIdx)
      powerRes$fdpKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 0)) / length(which(totOut$kernelAvg < 0.1))
      powerRes$powerKernel[sim_it] <- length(which(totOut$kernelAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
      powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
      powerRes$nRejKernel[sim_it] <-  length(which(totOut$kernelAvg < 0.1)) 
    } else {
      totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)
    }
  } else {totOut <- totOut %>% mutate(kernelLfdr = NA, kernelAvg=NA)}
  
  
  #   if (class(oldResKernel)[1] == "list") {
  #     totOut$kernelLfdr <- oldResKernel$lfdrVec
  #     totOut <- totOut %>% arrange(kernelLfdr) %>%
  #       mutate(kernelAvg = cummean(kernelLfdr)) %>%
  #       arrange(origIdx)
  #     rej_idx <- which(totOut$kernelAvg < 0.1)
  #     powerRes$fdpKernel[sim_it] <- if (length(rej_idx) > 0) sum(totOut$causal[rej_idx] == 0) / length(rej_idx) else NA
  #     powerRes$powerKernel[sim_it] <- if (sum(causalVec) > 0) sum(totOut$causal[rej_idx] == 1) / sum(causalVec) else NA
  #     powerRes$inconKernel[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldResKernel$lfdrVec))
  #     powerRes$nRejKernel[sim_it] <- length(rej_idx)
  #   }
  # }
  
  if (do7df) {
    oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                        dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=FALSE)
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
  
  #   if (class(oldRes7df)[1] == "list") {
  #     totOut$df7Lfdr <- oldRes7df$lfdrVec
  #     totOut <- totOut %>% arrange(df7Lfdr) %>%
  #       mutate(df7Avg = cummean(df7Lfdr)) %>%
  #       arrange(origIdx)
  #     rej_idx <- which(totOut$df7Avg < 0.1)
  #     powerRes$fdp7df[sim_it] <- if (length(rej_idx) > 0) sum(totOut$causal[rej_idx] == 0) / length(rej_idx) else NA
  #     powerRes$power7df[sim_it] <- if (sum(causalVec) > 0) sum(totOut$causal[rej_idx] == 1) / sum(causalVec) else NA
  #     powerRes$incon7df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes7df$lfdrVec))
  #     powerRes$nRej7df[sim_it] <- length(rej_idx)
  #   }
  # }
  
  if (do50df) {
    oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = allZ, kernel = FALSE, joint=FALSE, ind = TRUE,
                                         dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=FALSE)
    
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
  
  #   if (class(oldRes50df)[1] == "list") {
  #     totOut$df50Lfdr <- oldRes50df$lfdrVec
  #     totOut <- totOut %>% arrange(df50Lfdr) %>%
  #       mutate(df50Avg = cummean(df50Lfdr)) %>%
  #       arrange(origIdx)
  #     rej_idx <- which(totOut$df50Avg < 0.1)
  #     powerRes$fdp50df[sim_it] <- if (length(rej_idx) > 0) sum(totOut$causal[rej_idx] == 0) / length(rej_idx) else NA
  #     powerRes$power50df[sim_it] <- if (sum(causalVec) > 0) sum(totOut$causal[rej_idx] == 1) / sum(causalVec) else NA
  #     powerRes$incon50df[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = oldRes50df$lfdrVec))
  #     powerRes$nRej50df[sim_it] <- length(rej_idx)
  #   }
  # }
  
  # New method with correlation
  estCorAll <- cor(allZ)
  
  # # Update init lists for 3D
  # initPiListCor <- list(c(0.82), c(0.04), c(0.04), c(0.04), c(0.02), c(0.02), c(0.01), c(0.01))
  # # But simpler: use same structure as before but with 3D means
  # initMuListCor <- list(
  #   matrix(0, nrow = nDims, ncol = 1),
  #   matrix(c(3,0,0), nrow = nDims),
  #   matrix(c(0,3,0), nrow = nDims),
  #   matrix(c(0,0,3), nrow = nDims),
  #   matrix(c(3,3,0), nrow = nDims),
  #   matrix(c(3,0,3), nrow = nDims),
  #   matrix(c(0,3,3), nrow = nDims),
  #   matrix(c(8,8,8), nrow = nDims)
  # )
  
  initPiListCor <- list(c(0.82))
  for (i in 2:7) {initPiListCor[[i]] <- c(0.08 / 12, 0.08 / 12)}
  initPiListCor[[8]] <- c(0.1)
  # the csmGmm package will add the appropriate 0s to initMuList
  initMuListCor <- list(matrix(data=0, nrow=3, ncol=1))
  for (i in 2:7) {
    initMuListCor[[i]] <- cbind(rep(2, 3), rep(5, 3))
  }
  initMuListCor[[8]] <- matrix(data=c(8, 8, 8), nrow=3)
  
  
  corRes <- symm_fit_cor_EM_R1(t_value = t, testStats = as.matrix(allZ), corMat = estCorAll, 
                               initMuList = initMuListCor, initPiList = initPiListCor, eps = 1e-5, checkpoint = FALSE)
  
  # record
  totOut <- totOut %>% mutate(corLfdr = corRes$lfdrResults) %>%
    arrange(corLfdr) %>%
    mutate(corAvg = cummean(corLfdr)) %>%
    # don't forget to set it back to original index for further additions!!
    arrange(origIdx)
  powerRes$nRejNew[sim_it] <- length(which(totOut$corAvg < 0.1))
  powerRes$fdpNew[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 0)) /  length(which(totOut$corAvg < 0.1))
  powerRes$powerNew[sim_it] <- length(which(totOut$corAvg < 0.1 & totOut$causal == 1)) / sum(causalVec)
  powerRes$inconNew[sim_it] <- length(check_incongruous(zMatrix = allZ, lfdrVec = corRes$lfdrVec))
  
  
  
  cat('Done with sim', sim_it, '\n')
}

setwd(outputDir)
write.table(powerRes, outName, append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

