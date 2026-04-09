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
read_method_lfdr <- function(file_path, lfdr_col, avg_col, rej_col, n_ref, keep_only_col = NULL) {
  if (!file.exists(file_path)) {
    out <- data.frame(tmp = rep(NA_real_, n_ref))
    names(out) <- lfdr_col
  } else {
    raw <- fread(file_path)
    if (!is.null(keep_only_col) && keep_only_col %in% names(raw)) {
      out <- raw %>% transmute(!!lfdr_col := .data[[keep_only_col]])
    } else {
      out <- raw[, 1, drop = FALSE]
      names(out) <- lfdr_col
    }
  }
  
  out %>%
    mutate(idx = row_number()) %>%
    arrange(.data[[lfdr_col]]) %>%
    mutate(
      !!avg_col := cummean(.data[[lfdr_col]]),
      !!rej_col := as.integer(.data[[avg_col]] < qval)
    ) %>%
    arrange(idx) %>%
    select(-idx)
}

allResList <- list()
rejTab <- data.frame()

for (snp_it in 1:3) {
  tempRoot <- file.path(outputDir, paste0("med_analysis_aID", snp_it))
  tempStatsDat <- fread(paste0(tempRoot, "_dat.txt")) %>%
    mutate(SNP = snp_it)
  
  n_ref <- nrow(tempStatsDat)
  
  tempNewDat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_newlfdr.txt"),
    lfdr_col = "newLfdr",
    avg_col = "avgNew",
    rej_col = "rejNew",
    n_ref = n_ref,
    keep_only_col = "newLfdr"
  )
  
  tempKernelDat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_kernel.txt"),
    lfdr_col = "kernelLfdr",
    avg_col = "avgKernel",
    rej_col = "rejKernel",
    n_ref = n_ref
  )
  
  tempDf50Dat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_df50.txt"),
    lfdr_col = "Df50Lfdr",
    avg_col = "avgDf50",
    rej_col = "rejDf50",
    n_ref = n_ref
  )
  
  tempDf7Dat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_df7.txt"),
    lfdr_col = "Df7Lfdr",
    avg_col = "avgDf7",
    rej_col = "rejDf7",
    n_ref = n_ref
  )
  
  tempdeibDat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_deib.txt"),
    lfdr_col = "deibLfdr",
    avg_col = "avgdeib",
    rej_col = "rejdeib",
    n_ref = n_ref
  )
  
  tempmtestDat <- read_method_lfdr(
    file_path = paste0(tempRoot, "_mtest.txt"),
    lfdr_col = "mtestLfdr",
    avg_col = "avgmtest",
    rej_col = "rejmtest",
    n_ref = n_ref
  )
  
  fullDat <- bind_cols(
    tempStatsDat, tempNewDat, tempKernelDat, tempDf7Dat, tempDf50Dat, tempdeibDat, tempmtestDat
  )
  
  if (anyDuplicated(names(fullDat)) > 0) {
    dup_names <- names(fullDat)[duplicated(names(fullDat))]
    stop(paste("Duplicate columns in fullDat:", paste(unique(dup_names), collapse = ", ")))
  }
  
  allResList[[snp_it]] <- fullDat
  
  tempRej <- fullDat %>%
    filter(rejDf50 == 1 | rejDf7 == 1 | rejNew == 1 | rejKernel == 1 | rejdeib == 1 | rejmtest == 1)
  
  rejTab <- bind_rows(rejTab, tempRej)
  cat("done SNP:", snp_it, "\n")
}

allRej <- rejTab %>%
  mutate(numRej = rowSums(across(c(rejDf50, rejDf7, rejNew, rejKernel, rejdeib, rejmtest)), na.rm = TRUE)) %>%
  arrange(desc(numRej)) %>%
  slice(1:5) %>%
  select(Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej, SNP, rejNew, rejNew, rejNew,
         rejKernel, rejDf7, rejDf50, rejdeib, rejmtest)

allRej1 <- rejTab %>%
  mutate(numRej = rowSums(across(c(rejDf50, rejDf7, rejNew, rejKernel, rejdeib, rejmtest)), na.rm = TRUE)) %>%
  arrange(desc(numRej)) %>%
  select(Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej, SNP,
         avgNew, rejNew, avgKernel, rejKernel, avgDf7, rejDf7, avgDf50, rejDf50, avgdeib, rejdeib, avgmtest, rejmtest)

tab2DF <- data.frame(
  RS = c("rs55781567", "rs56113850", "rs7705526"),
  Gene = c("CHRNA5", "CYP2A6", "TERT"),
  BP = c(78857986, 41353107, 1285974),
  Chr = c(15, 19, 5),
  SNP = 1:3
)

# IMPORTANT: use left_join to avoid Gene.x / Gene.y collisions
mergedRej <- allRej %>%
  left_join(tab2DF %>% select(SNP, RS, Chr, BP), by = "SNP") %>%
  arrange(desc(numRej)) %>%
  select(RS, Chr, BP, Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej)

print(mergedRej)
