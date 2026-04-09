
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig6/mediation_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
# setwd("~/Downloads/csmGmm_sim_R3/Fig6")
here::i_am("Fig6/two_mediation_analysis.R")

# load libraries
library(mvtnorm)
library(usethis)
library(data.table)
library(dplyr)
library(Matrix)
library(expm)
library(magrittr)
library(devtools)
library(ks)
library(csmGmm)
library(here)
library(locfdr)
library(purrr)
library(DEIcompo)
library(mediation.test)

## basic functions
################################################################################
find_max_means_R1 <- function(muInfo) {
  
  # iterate, skip the first (0) and last (alternative)
  listLength <- length(muInfo)
  K <- nrow(muInfo[[1]])
  S1 <- c(1,2,4,5,6,8,9,10)+1
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
  # H_annot[which(
  #   t+1<=apply(abs(H_space), 1, sum) & 
  #     apply(abs(H_space), 1, sum) <= K &
  #     apply(H_space, 1, function(row) {
  #       (row[1]*row[2]) +(row[3]*row[4])
  #     })!=0
  # )] <- 1
  
  H_annot[which(
    t + 1 <= apply(abs(H_space), 1, sum) &
      apply(abs(H_space), 1, sum) <= K &
      apply(H_space, 1, function(row) {
        (row[1] != 0 && row[2] != 0) || (row[3] != 0 && row[4] != 0)
      })
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
  # if (sameDirAlt) {
  #   #null_rows <- which(apply(H_space, 1, sum) != ncol(summary_tab) & apply(H_space, 1, sum) != -ncol(summary_tab))
  #   null_rows <- which( apply(abs(H_space), 1, sum) <= t_value |
  #                        (apply(abs(H_space), 1, sum) > t_value &
  #                           apply(abs(H_space), 1, sum) <= ncol(summary_tab) &
  #                           apply(H_space, 1, function(row) {
  #                             row[1]*row[2] })==0 & apply(H_space, 1, function(row) { row[3]*row[4] })==0
  #                         ) )
  # } else {
  #   null_rows <- which(H_annot == 0)
  # }

  if (sameDirAlt) {
    nnz <- rowSums(abs(H_space))
    path12 <- (H_space[,1] != 0) & (H_space[,2] != 0)
    path34 <- (H_space[,3] != 0) & (H_space[,4] != 0)

    null_rows <- which(
      (nnz <= t_value) |
        ((nnz > t_value) & !path12 & !path34)
    )
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
  if (is.na(sum(neg_pmf_tab)) | is.na(sum(pos_pmf_tab))) {return(-1)}

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
  ## symmetric alternative
  # sum_12 <- apply(Hmat[, 1:K], 1, function(row) row[1] * row[2])
  # sum_34 <- apply(Hmat[, 1:K], 1, function(row) row[3] * row[4])
  # symAltVec <- ifelse((t+1 <= apply(abs(Hmat[, 1:K]), 1, sum)) & 
  #                       (apply(abs(Hmat[, 1:K]), 1, sum) <= K) &  
  #                       (sum_12 + sum_34 !=0) , 1, 0)
  #
  path_12 <- (Hmat[,1] != 0) & (Hmat[,2] != 0)
  path_34 <- (Hmat[,3] != 0) & (Hmat[,4] != 0)
  symAltVec <- ifelse(
    (t + 1 <= rowSums(abs(Hmat[, 1:K]))) &
      (rowSums(abs(Hmat[, 1:K])) <= K) &
      (path_12 | path_34),
    1, 0
  )
  
  
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
      as.data.frame() %>%
      #as.data.frame(.data) %>%
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

run_mtest_pair <- function(z1, z2, q = 0.1) {
  tt <- cbind(z1, z2)
  
  out <- suppressWarnings(
    tryCatch(
      mediation.test::mediation_test_minimax_BH(
        t = tt,
        alpha = q,
        BH = "statistics",
        sample_size = Inf
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(out) || is.null(out$decision)) {
    return(rep(NA_integer_, length(z1)))
  }
  
  as.integer(out$decision)
}


# helper: return minimax p-values (not decision) for one pair
run_mtest_pair_pval <- function(z1, z2, q = 0.1) {
  tt <- cbind(z1, z2)
  out <- suppressWarnings(
    tryCatch(
      mediation.test::mediation_test_minimax_BH(
        t = tt,
        alpha = q,
        BH = "pval",
        sample_size = Inf
      ),
      error = function(e) NULL
    )
  )
  if (is.null(out) || is.null(out$pvals)) return(rep(NA_real_, length(z1)))
  p <- out$pvals
  p[!is.finite(p)] <- 1
  p[is.na(p)] <- 1
  p
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
outputDir <- here::here("Fig6", "output")
outRoot <- paste0(outputDir, "/med_analysis_aID", aID)

# additional data needed
twasFname1 <- here::here("Data/mvp_scc__PM__Blood.csv")  #blood
#twasFname <- here::here("Data/mvp_scc__PM__Lung.csv")
twasFname <- here::here("Data/scc_lung_addchr1.csv")  #lung


# convergence of EM
oldEps <- 0.01
newEps <- 10^(-5)
t <- 1

# list of SNPs

tab2DF <- data.frame(RS = c("rs55781567", "rs56113850", "rs7705526"),
                     Gene = c("CHRNA5", "CYP2A6", "TERT"),
                     BP = c(78857986, 41353107, 1285974),
                     Chr = c(15, 19, 5)) %>%
  slice(aID:aID)

# read the lung twas results
twasRes <- read.csv(twasFname) %>%
  select(gene_name, zscore, pvalue) %>%
  set_colnames(c("Gene", "Z_twas", "p_twas")) %>%
  filter(!is.na(Z_twas))

# read the blood twas results
twasRes1 <- read.csv(twasFname1) %>%
  select(gene_name, zscore, pvalue) %>%
  set_colnames(c("Gene", "Z_twas1", "p_twas1")) %>%
  filter(!is.na(Z_twas1))

# loop through 15 SNPs
for (snp_it in 1:nrow(tab2DF)) {
  tempFname <- paste0(outputDir, "/", tab2DF$RS[snp_it], "_", tab2DF$Gene[snp_it], "_", 
                      tab2DF$Chr[snp_it], "_", tab2DF$BP[snp_it], ".txt")
  
  tempEqtl <- fread(tempFname) %>% 
    select(Gene, testStat, pval) %>%
    set_colnames(c("Gene", "Z_eqtl", "p_eqtl"))
  
  tempFname1 <- paste0(outputDir, "/", tab2DF$RS[snp_it], "blood", ".txt")
  
  tempEqtl1 <- fread(tempFname1) %>%
    mutate(pval1 = 1 - pchisq(Zscore^2, df = 1)) %>%
    select(GeneSymbol, Zscore, pval1) %>%
    set_colnames(c("Gene", "Z_eqtl1", "p_eqtl1"))
  
  # merge
  allDat <- tempEqtl %>% merge(twasRes, by="Gene") 
  allDat1 <- tempEqtl1 %>% merge(twasRes1, by="Gene") 
  #test statistics
  Combin_data <- allDat %>% merge(allDat1, by="Gene")
  
  cat("colnames", colnames(allDat1))
  
  allZ <-Combin_data %>%
    select(Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, p_twas, p_eqtl, p_twas1, p_eqtl1) %>%
    as.matrix(.)
  
  
  # adjustment to not get p-values of 0 needed for DEI-B and MTEST and our own methods
  for (col_it in 1:(ncol(allZ))/2) {
    tooBig <- which(allZ[, col_it] > 8.1)
    tooSmall <- which(allZ[, col_it] < -8.1)
    if (length(tooBig) > 0) {
      allZ[tooBig, col_it] <- 8.1
    }
    if (length(tooSmall) > 0) {
      allZ[tooSmall, col_it] <- -8.1
    }
  }
  
  for (col_it in ((ncol(allZ)/2)+1):ncol(allZ)) {
    tooBig <- which(allZ[, col_it] == 1)
    tooSmall <- which(allZ[, col_it] == 0)
    minVal <- min(allZ[which(allZ[, col_it] > 0), col_it])
    if (length(tooBig) > 0) {
      allZ[tooBig, col_it] <- 0.999
    }
    if (length(tooSmall) > 0) {
      allZ[tooSmall, col_it] <- minVal
    }
  }
  
  # save the test statistics
  write.table(Combin_data, paste0(outRoot, "_dat.txt"), append=F, quote=F, row.names=F, col.names=T)
  
  
  ## test stats / p-values
  testData <- allZ[, 1:4, drop = FALSE]
  pData    <- allZ[, 5:8, drop = FALSE]
  
  ## DEI-B (always run)
  p_dei_12 <- tryCatch(DEIcompo::DEI_B(pData[, 1], pData[, 2]), error = function(e) NULL)
  p_dei_34 <- tryCatch(DEIcompo::DEI_B(pData[, 3], pData[, 4]), error = function(e) NULL)
  
  if (is.null(p_dei_12) || is.null(p_dei_34) ||
      length(p_dei_12) != nrow(testData) || length(p_dei_34) != nrow(testData)) {
    deibLfdr <- rep(NA_real_, nrow(testData))
  } else {
    p_dei_12[!is.finite(p_dei_12)] <- 1
    p_dei_34[!is.finite(p_dei_34)] <- 1
    p_dei_12 <- pmin(pmax(p_dei_12, 0), 1)
    p_dei_34 <- pmin(pmax(p_dei_34, 0), 1)
    
    deibLfdr <- pmin(
      p.adjust(p_dei_12, method = "BH"),
      p.adjust(p_dei_34, method = "BH")
    )
  }
  
  write.table(deibLfdr, paste0(outRoot, "_deib.txt"),
              append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  ## MTEST (always run)
  p_m12 <- run_mtest_pair_pval(testData[, 1], testData[, 2], q = 0.1)
  p_m34 <- run_mtest_pair_pval(testData[, 3], testData[, 4], q = 0.1)
  
  if (all(is.na(p_m12)) && all(is.na(p_m34))) {
    mtestLfdr <- rep(NA_real_, nrow(testData))
  } else {
    p_m12[!is.finite(p_m12)] <- 1
    p_m34[!is.finite(p_m34)] <- 1
    p_m12 <- pmin(pmax(p_m12, 0), 1)
    p_m34 <- pmin(pmax(p_m34, 0), 1)
    
    mtestLfdr <- pmin(
      p.adjust(p_m12, method = "BH"),
      p.adjust(p_m34, method = "BH")
    )
  }
  
  write.table(mtestLfdr, paste0(outRoot, "_mtest.txt"),
              append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  
  # kernel
  # oldResKernel <- emp_bayes_framework(summary_tab = testData, sameDirAlt = FALSE, kernel = TRUE, joint=FALSE, ind = TRUE,
  #                                     dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)

  oldResKernel <- emp_bayes_framework_R1(t_value = t, summary_tab = testData, sameDirAlt=TRUE, kernel = TRUE, joint=FALSE, ind = TRUE,
                                         dfFit = 7, Hdist_epsilon=10^(-3), checkpoint=TRUE)


  if (class(oldResKernel)[1] != "list") {
    kernelLfdr <- rep(NA, nrow(testData))
  } else {
    kernelLfdr <- oldResKernel$lfdrVec
  }
  # save
  write.table(kernelLfdr, paste0(outRoot, "_kernel.txt"), append=F, quote=F, row.names=F, col.names=T)

  # 7 df
  # oldRes7df <- emp_bayes_framework(summary_tab = testData[,1:2], sameDirAlt = FALSE, kernel = FALSE, joint=FALSE, ind = TRUE,
  #                                  dfFit = 7, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  oldRes7df <- emp_bayes_framework_R1(t_value = t, summary_tab = testData, sameDirAlt=TRUE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                      dfFit = 7, Hdist_epsilon=10^(-3), checkpoint=TRUE)

  if (class(oldRes7df)[1] != "list") {
    df7Lfdr <- rep(NA, nrow(testData))
  } else {
    df7Lfdr <- oldRes7df$lfdrVec
  }
  # save
  write.table(df7Lfdr, paste0(outRoot, "_df7.txt"), append=F, quote=F, row.names=F, col.names=T)

  # 50 df
  # oldRes50df <- emp_bayes_framework(summary_tab = testData, sameDirAlt = FALSE, kernel = FALSE, joint=FALSE, ind = TRUE,
  #                                   dfFit = 50, Hdist_epsilon=10^(-2), checkpoint=TRUE)
  oldRes50df <- emp_bayes_framework_R1(t_value = t, summary_tab = testData, sameDirAlt=TRUE, kernel = FALSE, joint=FALSE, ind = TRUE,
                                       dfFit = 50, Hdist_epsilon=10^(-3), checkpoint=TRUE)

  if (class(oldRes50df)[1] != "list") {
    df50Lfdr <- rep(NA, nrow(testData))
  } else {
    df50Lfdr <- oldRes50df$lfdrVec
  }
  # save
  write.table(df50Lfdr, paste0(outRoot, "_df50.txt"), append=F, quote=F, row.names=F, col.names=T)

  
  
  # new method
  initPiList <- list(c(0.82))
  for (i in 2:(2^4)) {initPiList[[i]] <- 0.18 / (2^4 - 1)}
  # the symm_fit_ind.R code will add the appropriate 0s to initMuList
  initMuList <- list(matrix(data=rep(0, 4), nrow=4, ncol=1))
  for (i in 2:(2^4)) {
    initMuList[[i]] <- matrix(data=rep(3, 4), nrow=4, ncol=1)
  }
  
  # nDims <- ncol(testData)
  # initPiList <- list(c(0.82))
  # for (i in 2:(2^nDims-1)) {initPiList[[i]] <- c(0.08 / 28, 0.08 / 28)}
  # initPiList[[2^nDims]] <- c(0.1)
  # # the symm_fit_ind.R code will add the appropriate 0s to initMuList
  # initMuList <- list(matrix(data=0, nrow=nDims, ncol=1))
  # for (i in 2:(2^nDims-1)) {
  #   initMuList[[i]] <- cbind(rep(2, nDims), rep(5, nDims))
  # }
  # initMuList[[2^nDims]] <- matrix(data=c(8, 8, 8, 8), nrow=nDims)
  
  
  #newRes <- symm_fit_ind_EM(testStats = testDat[, 1:2], sameDirAlt = FALSE, initMuList = initMuList, initPiList = initPiList, eps=newEps)
  
  #corRes <- symm_fit_cor_EM(t_value = t, testStats = testData, corMat = estCorAll, sameDirAlt = TRUE, initMuList = initMuListCor, initPiList = initPiListCor, eps=newEps)
  newRes <- symm_fit_ind_EM_R1(t_value = t, testStats = testData, initMuList = initMuList,
                               initPiList = initPiList, sameDirAlt=TRUE, eps=10^(-4))
  
  # save
  newOut <- Combin_data %>% mutate(newLfdr = newRes$lfdrResults)
  write.table(newOut, paste0(outRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
  write.table(do.call(cbind, newRes$muInfo), paste0(outRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
  write.table(do.call(cbind, newRes$piInfo), paste0(outRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
}


