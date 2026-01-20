# Process raw analysis of UKB data

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig7/summarize_ukb_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
here::i_am("Fig7/summarize_ukb_analysis_R1.R")

# load libraries
library(data.table)
library(dplyr)
library(devtools)
library(ks)
library(csmGmm)
library(here)

# record input - controls seed, parameters, etc.
args <- commandArgs(trailingOnly=TRUE)
aID <- as.numeric(args[1])
Snum <- as.numeric(args[2])

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Fig7", "output")
dataDir <- here::here("Data")
fnameOut <- paste0(outputDir, "/processed_ukb_data_S", Snum, ".txt")
rejectFnameRoot <- paste0(outputDir, "/reject_bmi_with_overall_neg5_reject_S_", Snum)

# nominal fdr
if (Snum == 1) {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.01
  fdrLimitKerneli <- 0.01
  fdrLimit7i <- 0.01
  fdrLimit50i <- 0.01
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.01
  fdrLimitDACTr <- 0.01
  fdrLimitKernelr <- 0.01
  fdrLimit7r <- 0.01
  fdrLimit50r <- 0.01
} else {
  fdrLimitHDMTi <- 0.1
  fdrLimitDACTi <- 0.1
  fdrLimitKerneli <- 0.1
  fdrLimit7i <- 0.1
  fdrLimit50i <- 0.1
  fdrLimitNew <- 0.1
  fdrLimitHDMTr <- 0.1
  fdrLimitDACTr <- 0.1
  fdrLimitKernelr <- 0.1
  fdrLimit7r <- 0.1
  fdrLimit50r <- 0.1
}


# how the raw output files are named
fnameRoot <- paste0(outputDir, "/Fig4_data_aID", 1:2)
fnameKernel <- paste0(fnameRoot, "_kernel.txt")
fname7 <- paste0(fnameRoot, "_df7.txt")
fname50 <- paste0(fnameRoot, "_df50.txt")
fnameNew <- c(paste0(fnameRoot, "_newlfdr.txt"))

selections <- list()
selections[[1]] <- c("Zoverall", "Zcad", "Zbmi")
selections[[2]] <- c("zILCCO", "zUKB", "zMVP")

# results
allResults <- c()
for (file_it in 1:2) {
  
  if (file_it == 1) {
    cleanZ <- fread(here::here(dataDir, "bmi_with_overall.txt"))
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7 <- fdrLimit7i 
  } else {
    cleanZ <- fread(here::here(dataDir, "lc_overall_three.txt"))
    fdrLimitKernel <- fdrLimitKerneli
    fdrLimit50 <- fdrLimit50i
    fdrLimit7 <- fdrLimit7i  
  }
  
  
  # hold temporary results
  tempRes <- data.frame(Method=c("New", "Kernel", "df50", "df7"), numReject=NA)
  
  # new 
  tempNew <- fread(fnameNew[file_it], header=T, data.table=F)
  tempDat <- cleanZ %>% select(all_of(selections[[file_it]]), chrpos) %>%
    mutate(origIdx = 1:nrow(.)) %>%   
    mutate(newLfdr = tempNew$x) %>%
    arrange(newLfdr) %>%
    mutate(cumNew = cummean(newLfdr)) %>%
    mutate(rejNew = ifelse(cumNew < fdrLimitNew, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[1] <- sum(tempDat$rejNew)
  
  
  
  # kernel
  tempKernel <- fread(fnameKernel[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(kernelLfdr = tempKernel$x) %>%
    arrange(kernelLfdr) %>%
    mutate(cumKernel = cummean(kernelLfdr)) %>%
    mutate(rejKernel = ifelse(cumKernel < fdrLimitKernel, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[2] <- sum(tempDat$rejKernel)
  
  # df50
  tempdf50 <- fread(fname50[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df50Lfdr = tempdf50$x) %>%
    arrange(df50Lfdr) %>%
    mutate(cumdf50 = cummean(df50Lfdr)) %>%
    mutate(rejdf50 = ifelse(cumdf50 < fdrLimit50, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[3] <- sum(tempDat$rejdf50)
  
  # df7
  tempdf7 <- fread(fname7[file_it], header=T, data.table=F)
  tempDat <- tempDat %>% mutate(df7Lfdr = tempdf7$x) %>%
    arrange(df7Lfdr) %>%
    mutate(cumdf7 = cummean(df7Lfdr)) %>%
    mutate(rejdf7 = ifelse(cumdf7 < fdrLimit7, 1, 0)) %>%
    arrange(origIdx)
  tempRes$numReject[4] <- sum(tempDat$rejdf7)
  
  
  # any rejection
  tempDat <- tempDat %>% mutate(rejAny = ifelse(rejdf7 == 1 | rejdf50 == 1 | rejKernel == 1 | rejNew == 1, 1, 0))
  
  rejectDat <- tempDat %>% filter(rejAny == 1) 
  
  # allResults
  allResults <- rbind(allResults, tempRes %>% mutate(aID = file_it))
  
  # save rejections
  write.table(rejectDat, paste0(rejectFnameRoot, "_aID", file_it, ".txt"), append=F, quote=F, row.names=F, col.names=T, sep='\t')
  
  cat(file_it)
}

# save
write.table(allResults, fnameOut, append=F, quote=F, row.names=F, col.names=T, sep='\t')






