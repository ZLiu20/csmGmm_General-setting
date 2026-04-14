# Make Figure 6 and Tables 1-2 

# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig4/plot_data_analysis.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.
setwd("~/Downloads/csmGmm_sim_R3/Fig6")
here::i_am("Fig6/plot_data_analysis.R")

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
outputDir <- here::here("Fig6", "output")
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

ggsave(paste0(outputDir, "/Fig_Pleiotropy.pdf"), width=14, height=18)


#----------------------------------------------------------------------------#
# Table 2

adjAnal  <- fread(here::here(outputDir, "processed_ukb_data_S1.txt"))
origAnal <- fread(here::here(outputDir, "processed_ukb_data_S2.txt"))

adaFilter <- fread(here::here(outputDir, "Fig4_data_aID1_adaFilter.txt"))
adjqch    <- fread(here::here(outputDir, "Fig4_data_aID1_qch01.txt"))
origqch   <- fread(here::here(outputDir, "Fig4_data_aID1_qch10.txt"))

# -----------------------------
# Harmonize method names
# -----------------------------
map_method <- function(x) {
  fcase(
    x == "New",  "csmGmm",
    x == "df50", "locfdr50",
    x == "df7",  "locfdr7",
    default = x
  )
}
origAnal[, Method := map_method(Method)]
adjAnal[,  Method := map_method(Method)]

# -----------------------------
# Compute extra method counts (aID=1 only)
# -----------------------------
num_rejadaFilter <- sum(unlist(adaFilter) < 0.1,  na.rm = TRUE)
num_rejQCH_adj   <- sum(unlist(adjqch)    < 0.01, na.rm = TRUE)
num_rejQCH_orig  <- sum(unlist(origqch)   < 0.1,  na.rm = TRUE)

orig_extra <- data.table(
  Method = c("qch_copula", "adaFilter"),
  numReject = c(num_rejQCH_orig, num_rejadaFilter),
  aID = 1L
)
adj_extra <- data.table(
  Method = c("qch_copula", "adaFilter"),
  numReject = c(num_rejQCH_adj, num_rejadaFilter),
  aID = 1L
)

# -----------------------------
# Combine tables
# -----------------------------
orig_all <- rbindlist(list(origAnal, orig_extra), use.names = TRUE, fill = TRUE)
adj_all  <- rbindlist(list(adjAnal,  adj_extra),  use.names = TRUE, fill = TRUE)

orig_all[, Source := "Orig"]
adj_all[,  Source := "Adj"]

long <- rbindlist(list(orig_all, adj_all), use.names = TRUE, fill = TRUE)[
  !is.na(Method) & !is.na(aID),
  .(Method = as.character(Method), aID = as.integer(aID), Source, numReject = as.integer(numReject))
]

# aID meaning
long[, Analysis := fcase(
  aID == 1L, "Pleiotropy",
  aID == 2L, "Replication",
  default = paste0("aID_", aID)
)]

# -----------------------------
# Build Table 2
# -----------------------------
tab2final <- dcast(
  long[aID %in% c(1L, 2L)],
  Method ~ Source + Analysis,
  value.var = "numReject",
  fun.aggregate = function(x) x[1]
)

# ensure expected columns exist
need_cols <- c(
  "Orig_Pleiotropy", "Adj_Pleiotropy",
  "Orig_Replication", "Adj_Replication"
)
for (nm in need_cols) {
  if (!nm %in% names(tab2final)) tab2final[, (nm) := NA_integer_]
}

setcolorder(
  tab2final,
  c("Method", "Orig_Pleiotropy", "Adj_Pleiotropy", "Orig_Replication", "Adj_Replication")
)

# optional display order
method_order <- c("csmGmm", "Kernel", "locfdr50", "locfdr7", "qch_copula", "adaFilter")
tab2final[, ord := fifelse(Method %in% method_order, match(Method, method_order), 999L)]
setorder(tab2final, ord, Method)
tab2final[, ord := NULL]

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

# print(mergedRej)

mergedRej1 <- allRej1 %>%
  left_join(tab2DF %>% select(SNP, RS), by = "SNP") %>%
  arrange(desc(numRej)) %>%
  select(RS, Gene, Z_eqtl, Z_twas, Z_eqtl1, Z_twas1, numRej,
         avgNew, avgKernel, avgDf7,avgdeib, avgmtest)

mergedRej1 <- mergedRej1 %>% tidyr::drop_na()


# ----------------------------------------------------------------------------------------
# Create all possible pairs of rows
result <- mergedRej1 %>%
  # Add row numbers for tracking
  mutate(row_id = row_number()) %>%
  # Create all pairs using cross join
  cross_join(mergedRej1 %>% mutate(row_id = row_number()), 
             suffix = c("_1", "_2")) %>%
  # Ensure we only compare each pair once (avoid duplicates and self-comparisons)
  filter(row_id_1 < row_id_2) %>%
  # Apply your filtering conditions
  filter(
    # Same RS values
    RS_1 == RS_2 &
      
      # Same signs for all Z values between the two rows
      sign(Z_eqtl_1) == sign(Z_eqtl_2) &
      sign(Z_eqtl1_1) == sign(Z_eqtl1_2) &
      sign(Z_twas_1) == sign(Z_twas_2) &
      sign(Z_twas1_1) == sign(Z_twas1_2) &
      
      # Either all abs values from row 1 are larger OR all are smaller than row 2
      (
        (abs(Z_eqtl_1) > abs(Z_eqtl_2) & 
           abs(Z_eqtl1_1) > abs(Z_eqtl1_2) & 
           abs(Z_twas_1) > abs(Z_twas_2) & 
           abs(Z_twas1_1) > abs(Z_twas1_2) ) |
          (abs(Z_eqtl_1) < abs(Z_eqtl_2) & 
             abs(Z_eqtl1_1) < abs(Z_eqtl1_2) & 
             abs(Z_twas_1) < abs(Z_twas_2) & 
             abs(Z_twas1_1) < abs(Z_twas1_2) ) 
      ) &
      
      (
        (avgKernel_1 < 0.1 & avgKernel_2 > 0.1) | (avgKernel_1 > 0.1 & avgKernel_2 < 0.1) |
          (avgDf7_1 < 0.1 & avgDf7_2 > 0.1) | (avgDf7_1 > 0.1 & avgDf7_2 < 0.1)
      ) &
      
      (
        (avgKernel_1 > avgKernel_2 & avgDf7_1 > avgDf7_2 ) |
          (avgKernel_1 < avgKernel_2 & avgDf7_1 < avgDf7_2 )
      ) &
      
      # Keep the numRej condition if needed (applied to both rows)
      numRej_1 < 3 & numRej_2 < 3
  ) %>%
  # Reshape to show pairs in two separate rows
  pivot_longer(
    cols = -c(row_id_1, row_id_2),
    names_to = c(".value", "pair"),
    names_pattern = "(.+)_(.)"
  ) %>%
  # Round Z_eqtl to 3 decimal places
  mutate(
    Z_eqtl = round(Z_eqtl, 2),
    Z_eqtl1 = round(Z_eqtl1, 2),
    Z_twas = round(Z_twas, 2),
    Z_twas1 = round(Z_twas1, 2),
    # avgNew = round(avgNew, 6),
    # avgKernel = round(avgKernel, 3),
    # avgDf7 = round(avgDf7, 3)
    
  ) %>%
  arrange(row_id_1, row_id_2) %>%
  distinct() %>%
  select(-c(1:3)) 


# Display the final result
print(result)

# Print the selected genes
print(mergedRej)

