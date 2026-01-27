# Make Figure 7 and Table 1 

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("~/Downloads/csmGmm_sim_R3/Fig7")
here::i_am("Fig7/plot_data_analysis.R")

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)
library(devtools)
library(here)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# set output directory 
outputDir <- here::here("Fig7", "output")
dataDir <- here::here("Data")

# for colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))

  # add true positions
  truePos <- rep(NA, nrow(plotRes))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotRes %>% filter(Chr == tempChr)
    truePos[counter:(counter + nrow(tempDat) - 1)] <- rep(sum(chrCounts[1:tempChr]), nrow(tempDat)) + tempDat$BP
    counter <- counter + nrow(tempDat)
  }

  # plot
  xBreaks <- cumsum(chrCounts[-1])
  xBreaksLabs <- 1:22
  #xBreaksLabs[c(9, 11, 13, 15, 16, 17, 19, 20, 21)] <- ""

  print(xBreaks)
  print(xBreaksLabs)
  
  plotDat <- plotRes %>% mutate(truePos = truePos)

  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))


  return(returnPlot)
}

# add position information to data
s1 <- fread(here::here(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID1.txt"))
s2 <- fread(here::here(outputDir, "reject_bmi_with_overall_neg5_reject_S_1_aID2.txt"))

s1new <- s1 %>% filter(rejNew == 1) %>%
  #mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))

# s1new <- s1 %>% filter(rejNew == 1) %>%
#   mutate(chars = nchar(chrpos)) %>%
#   mutate(colonPos = unlist(gregexpr(":", chrpos))) %>%
#   mutate(Chr = substr(chrpos, 1, colonPos - 1),
#          BP = substr(chrpos, colonPos + 1, chars)) %>%
#   mutate(Chr = as.numeric(Chr), BP = as.numeric(BP))  


s2new <- s2 %>% filter(rejNew == 1) %>%
  #mutate(chars = nchar(chrpos)) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1))  %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars))


# for plotting axes
# allZ <- fread(here::here(dataDir, "bmi_with_overall.txt")) %>%
allZ <- fread(here::here(dataDir, "lc_overall_three.txt")) %>%
  mutate(chars = nchar(chrpos)) %>%
  mutate(colonPos = gregexpr(":", chrpos)) %>%
  mutate(colonPos = as.numeric(colonPos)) %>%
  mutate(Chr = substr(chrpos, 1, colonPos - 1)) %>%
  mutate(BP = substr(chrpos, colonPos + 1, chars)) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP))
chrCounts <- rep(0, 23)
for (chr_it in 1:22) {
  tempDat <- allZ %>% filter(Chr == chr_it)
  maxPos <- max(tempDat$BP)
  chrCounts[chr_it + 1] <- maxPos
}

# data for manhattan plot - 2 way pleiotropy
# data for manhattan plot - Pleiotropy & replication
manDataRep <- rbind(s1new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "Pleiotropy"),
                    s2new %>% select(Chr, BP, chrpos, newLfdr) %>% mutate(pheno = "Replication")) %>%
  as.data.frame(.) %>%
  mutate(Chr = as.numeric(Chr), BP = as.numeric(BP)) %>%
  # overlap with LC, CAD pleiotropy
  mutate(Pleio = ifelse(chrpos %in% s1new$chrpos, 1, 0)) %>%
  mutate(Rep = ifelse(chrpos %in% s2new$chrpos, 1, 0)) %>%
  mutate(cat = ifelse(Pleio == 1, "Pleiotropy (CAD,BMI,LC)","Replication (ILCCO, UKB, MVP)")) %>%
  arrange(pheno) %>%
  # distinct() keeps the first one, so we arrange first - want to show Pleio lfdr if in pleio
  distinct(., chrpos, .keep_all = TRUE)

manPlotRep <- plotManhattan(plotRes = manDataRep, chrCounts,
                            colValues=c(gg_color_hue(3)[3], "darkorange"), shapeValues=c(17, 18),
                            ylimits=c(0, 6.5), legName="Pleiotropy & Replication")
manPlotTwo

ggsave(paste0(outputDir, "/Fig_Pleiotropy.pdf"), width=18, height=8)


#----------------------------------------------------------------------------#
# Table 2

# read summary of pleiotropy analysis
adjAnal <- fread(here::here(outputDir, "processed_ukb_data_S1.txt"))
origAnal <- fread(here::here(outputDir, "processed_ukb_data_S2.txt"))

tab2 <- origAnal %>% filter(aID == 1) %>% select(Method, numReject) %>%
  mutate(a1 = adjAnal %>% filter(aID == 1) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a1new = origAnal %>% filter(aID == 2) %>% select(numReject) %>% unlist(.)) %>%
  mutate(a2new = adjAnal %>% filter(aID == 2) %>% select(numReject) %>% unlist(.)) 

tab2final <- tab2 %>%
  set_colnames(c("Method", "Orig", "Adj", "Orig ", "Adj ")) %>%
  mutate(Method = c("csmGmm", "Kernel", "locfdr50", "locfdr7"))
tab2final


#----------------------------------------------------------------------------------------#


# Table 1
qval <- 0.1

# initDat <- fread(here::here(outputDir, "med_analysis_aID1_newlfdr.txt")) %>%
#   select(Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1)

# hold results
allResList <- list()
rejTab <- c()
for(snp_it in setdiff(1:15,c(1,4,12,15))) {
  
  tempRoot <- paste0(outputDir, "/med_analysis_aID", snp_it)
  tempStatsDat <- fread(paste0(tempRoot, "_dat.txt")) %>%
    mutate(SNP = snp_it)
  
  
  if (file.exists(paste0(tempRoot, "_newlfdr.txt"))) { 
    tempNewDat <- fread(paste0(tempRoot, "_newlfdr.txt")) %>%
      mutate(idx = 1:nrow(.))  %>%
      arrange(newLfdr) %>%
      mutate(avgNew = cummean(newLfdr)) %>%
      mutate(rejNew = ifelse(avgNew < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx, -Gene, -Z_eqtl, -Z_twas, -Z_eqtl1, -Z_twas1)
  } else {tempNewDat <- data.frame(newLfdr=NA, rejNew=rep(0, nrow(tempStatsDat)), avgNew=NA)}
  
  if (file.exists(paste0(tempRoot, "_kernel.txt"))) {
    tempKernelDat <- fread(paste0(tempRoot, "_kernel.txt")) %>% set_colnames(c("kernelLfdr")) %>%
      mutate(idx = 1:nrow(.)) %>% 
      arrange(kernelLfdr) %>%
      mutate(avgKernel = cummean(kernelLfdr)) %>%
      mutate(rejKernel = ifelse(avgKernel < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx) 
  } else {tempKernelDat <- data.frame(kernelLfdr=NA, rejKernel=rep(0, nrow(tempStatsDat)), avgKernel=NA)}
  
  if (file.exists(paste0(tempRoot, "df50.txt"))) { 
    tempDf50Dat <- fread(paste0(tempRoot, "df50.txt")) %>% set_colnames(c("Df50Lfdr")) %>%
      mutate(idx = 1:nrow(.)) %>%      
      arrange(Df50Lfdr) %>%
      mutate(avgDf50 = cummean(Df50Lfdr)) %>%
      mutate(rejDf50 = ifelse(avgDf50 < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx) 
  } else {tempDf50Dat <- data.frame(Df50Lfdr=NA, rejDf50=rep(0, nrow(tempStatsDat)), avgDf50=NA)}
  
  if (file.exists(paste0(tempRoot, "_df7.txt"))) {
    tempDf7Dat <- fread(paste0(tempRoot, "_df7.txt")) %>% set_colnames(c("Df7Lfdr")) %>%
      mutate(idx = 1:nrow(.)) %>%      
      arrange(Df7Lfdr) %>%
      mutate(avgDf7 = cummean(Df7Lfdr)) %>%
      mutate(rejDf7 = ifelse(avgDf7 < qval, 1, 0)) %>%
      arrange(idx) %>%
      select(-idx)   
  } else {tempDf7Dat <- data.frame(Df7Lfdr=NA, rejDf7=rep(0, nrow(tempStatsDat)), avgDf7=NA)} 
  
  
  # put in list
  fullDat <- cbind(tempStatsDat, tempNewDat, tempKernelDat, tempDf7Dat)  # df50 method can't work
  allResList[[snp_it]] <- fullDat
  
  # rejections
  tempRej <- fullDat %>% filter( rejDf7 == 1 | rejNew == 1 | rejKernel == 1)
  rejTab <- rbind(rejTab, tempRej)
  
  cat(snp_it)
}

# look at rejections
allRej <- rejTab %>% mutate(numRej = rejDf7 + rejNew + rejKernel, na.rm=TRUE) %>%
  arrange(desc(numRej)) %>%
  slice(1:5) %>%
  select(Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej, SNP)

# merge with information
tab2DF <- data.frame(RS = c("rs71658797", "rs6920364", "rs11780471",
                            "rs55781567", "rs56113850", "rs13080835",
                            "rs7705526", "rs4236709", "rs885518", "rs11591710",
                            "rs1056562", "rs77468143", "rs41309931", "rs116822326",
                            "rs7953330"),
                     Locus = c("FUBP1", "RNASET2", "CHRNA2",
                               "CHRNA5", "CYP2A6", "TP63", "TERT", "NRG1", "CDKN2A",
                               "OBFC1", "AMICA1", "SECISBP2L", "RTEL1", "HCP5", "RAD52"),
                     BP = c(77967507, 167376466, 27344719, 78857986,
                            41353107, 189357199, 1285974, 32410110, 21830157, 105687632,
                            118125625, 49376624, 62326579, 31434111, 998819),
                     Chr = c(1, 6, 8, 15, 19, 3, 5, 8, 9, 10, 11, 15, 20, 6, 12), SNP=1:15)

mergedRej <- merge(allRej, tab2DF, by="SNP") %>%
  select(RS, Chr, BP, Locus, Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej)

mergedRej
