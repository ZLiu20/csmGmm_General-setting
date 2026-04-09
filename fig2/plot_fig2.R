# Collect results and plot Figure 2
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig2/plot_fig2.R") or set the path after the -cwd flag
# in the .lsf file, and then run again

setwd("~/Downloads/csmGmm_sim_R3/Fig2")
here::i_am("Fig2/plot_fig2.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(csmGmm)
library(here)


summarize_raw_R1 <- function(fullDat, full = FALSE, cor = FALSE, maxP = FALSE, FDP2 = FALSE) {
  safe_mean <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  }
  
  out_list <- list()
  effSizes <- sort(unique(fullDat$minEff1))
  
  for (i in seq_along(effSizes)) {
    tempEff <- effSizes[i]
    
    tempDat <- fullDat %>%
      dplyr::filter(minEff1 == tempEff) %>%
      as.data.frame() %>%
      dplyr::mutate(
        fdpDACT   = ifelse(nRejDACT   == 0, 0, fdpDACT),
        fdpHDMT   = ifelse(nRejHDMT   == 0, 0, fdpHDMT),
        fdpKernel = ifelse(nRejKernel == 0, 0, fdpKernel),
        fdp7df    = ifelse(nRej7df    == 0, 0, fdp7df),
        fdp50df   = ifelse(nRej50df   == 0, 0, fdp50df),
        fdpNew    = ifelse(nRejNew    == 0, 0, fdpNew),
        fdpDEIB   = ifelse(nRejDEIB   == 0, 0, fdpDEIB),
        fdpMTEST  = ifelse(nRejMTEST  == 0, 0, fdpMTEST)
      )
    
    if (cor) {
      tempDat <- tempDat %>%
        dplyr::mutate(
          fdpCor = ifelse(nRejCor == 0, 0, fdpCor),
          fdpNew = fdpCor, powerNew = powerCor, nRejNew = nRejCor
        )
    }
    
    if (full) {
      tempDat <- tempDat %>%
        dplyr::mutate(
          fdpFull = ifelse(nRejFull == 0, 0, fdpFull),
          fdpDACT = fdpFull, powerDACT = powerFull, nRejDACT = nRejFull
        )
    }
    
    if (maxP) {
      methods <- c("MaxP", "DACT", "HDMT", "Kernel", "df7", "df50", "New", "DEIB", "MTEST")
      
      summaryOut <- data.frame(
        minEff1 = tempDat$minEff1[1],
        Method = methods,
        stringsAsFactors = FALSE
      )
      
      summaryOut$nRej <- c(
        safe_mean(tempDat$nRejMaxp), safe_mean(tempDat$nRejDACT), safe_mean(tempDat$nRejHDMT),
        safe_mean(tempDat$nRejKernel), safe_mean(tempDat$nRej7df), safe_mean(tempDat$nRej50df),
        safe_mean(tempDat$nRejNew), safe_mean(tempDat$nRejDEIB), safe_mean(tempDat$nRejMTEST)
      )
      
      summaryOut$Power <- c(
        safe_mean(tempDat$powerMaxp), safe_mean(tempDat$powerDACT), safe_mean(tempDat$powerHDMT),
        safe_mean(tempDat$powerKernel), safe_mean(tempDat$power7df), safe_mean(tempDat$power50df),
        safe_mean(tempDat$powerNew), safe_mean(tempDat$powerDEIB), safe_mean(tempDat$powerMTEST)
      )
      
      summaryOut$FDP <- c(
        safe_mean(tempDat$fdpMaxp), safe_mean(tempDat$fdpDACT), safe_mean(tempDat$fdpHDMT),
        safe_mean(tempDat$fdpKernel), safe_mean(tempDat$fdp7df), safe_mean(tempDat$fdp50df),
        safe_mean(tempDat$fdpNew), safe_mean(tempDat$fdpDEIB), safe_mean(tempDat$fdpMTEST)
      )
      
      summaryOut$Incongruous <- c(
        NA_real_, NA_real_, NA_real_,
        safe_mean(tempDat$inconKernel), safe_mean(tempDat$incon7df),
        safe_mean(tempDat$incon50df), safe_mean(tempDat$inconNew),
        safe_mean(tempDat$inconDEIB), safe_mean(tempDat$inconMTEST)
      )
      
      summaryOut$numNA <- c(
        sum(is.na(tempDat$powerMaxp)), sum(is.na(tempDat$powerDACT)), sum(is.na(tempDat$powerHDMT)),
        sum(is.na(tempDat$powerKernel)), sum(is.na(tempDat$power7df)), sum(is.na(tempDat$power50df)),
        sum(is.na(tempDat$powerNew)), sum(is.na(tempDat$powerDEIB)), sum(is.na(tempDat$powerMTEST))
      )
      
      if (FDP2) {
        summaryOut$sig2Pow <- c(
          NA_real_, safe_mean(tempDat$sig2PowDACT), safe_mean(tempDat$sig2PowHDMT),
          safe_mean(tempDat$sig2PowKernel), safe_mean(tempDat$sig2Pow7df),
          safe_mean(tempDat$sig2Pow50df), safe_mean(tempDat$sig2PowNew),
          safe_mean(tempDat$sig2PowDEIB), safe_mean(tempDat$sig2PowMTEST)
        )
        
        summaryOut$sig2FDb.P <- c(
          NA_real_, safe_mean(tempDat$sig2fdpDACT), safe_mean(tempDat$sig2fdpHDMT),
          safe_mean(tempDat$sig2fdpKernel), safe_mean(tempDat$sig2fdp7df),
          safe_mean(tempDat$sig2fdp50df), safe_mean(tempDat$sig2fdpNew),
          safe_mean(tempDat$sig2fdpDEIB), safe_mean(tempDat$sig2fdpMTEST)
        )
      }
      
    } else {
      methods <- c("DACT", "HDMT", "Kernel", "df7", "df50", "New", "DEIB", "MTEST")
      
      summaryOut <- data.frame(
        minEff1 = tempDat$minEff1[1],
        Method = methods,
        stringsAsFactors = FALSE
      )
      
      summaryOut$nRej <- c(
        safe_mean(tempDat$nRejDACT), safe_mean(tempDat$nRejHDMT), safe_mean(tempDat$nRejKernel),
        safe_mean(tempDat$nRej7df), safe_mean(tempDat$nRej50df), safe_mean(tempDat$nRejNew),
        safe_mean(tempDat$nRejDEIB), safe_mean(tempDat$nRejMTEST)
      )
      
      summaryOut$Power <- c(
        safe_mean(tempDat$powerDACT), safe_mean(tempDat$powerHDMT), safe_mean(tempDat$powerKernel),
        safe_mean(tempDat$power7df), safe_mean(tempDat$power50df), safe_mean(tempDat$powerNew),
        safe_mean(tempDat$powerDEIB), safe_mean(tempDat$powerMTEST)
      )
      
      summaryOut$FDP <- c(
        safe_mean(tempDat$fdpDACT), safe_mean(tempDat$fdpHDMT), safe_mean(tempDat$fdpKernel),
        safe_mean(tempDat$fdp7df), safe_mean(tempDat$fdp50df), safe_mean(tempDat$fdpNew),
        safe_mean(tempDat$fdpDEIB), safe_mean(tempDat$fdpMTEST)
      )
      
      summaryOut$Incongruous <- c(
        NA_real_, NA_real_, safe_mean(tempDat$inconKernel), safe_mean(tempDat$incon7df),
        safe_mean(tempDat$incon50df), safe_mean(tempDat$inconNew), 
        safe_mean(tempDat$inconDEIB), safe_mean(tempDat$inconMTEST)
      )
      
      summaryOut$numNA <- c(
        sum(is.na(tempDat$powerDACT)), sum(is.na(tempDat$powerHDMT)), sum(is.na(tempDat$powerKernel)),
        sum(is.na(tempDat$power7df)), sum(is.na(tempDat$power50df)), sum(is.na(tempDat$powerNew)),
        sum(is.na(tempDat$powerDEIB)), sum(is.na(tempDat$powerMTEST))
      )
      
      if (FDP2) {
        summaryOut$sig2Pow <- c(
          safe_mean(tempDat$sig2PowDACT), safe_mean(tempDat$sig2PowHDMT), safe_mean(tempDat$sig2PowKernel),
          safe_mean(tempDat$sig2Pow7df), safe_mean(tempDat$sig2Pow50df), safe_mean(tempDat$sig2PowNew),
          safe_mean(tempDat$sig2PowDEIB), safe_mean(tempDat$sig2PowMTEST)
        )
        
        summaryOut$sig2FDb.P <- c(
          safe_mean(tempDat$sig2fdpDACT), safe_mean(tempDat$sig2fdpHDMT), safe_mean(tempDat$sig2fdpKernel),
          safe_mean(tempDat$sig2fdp7df), safe_mean(tempDat$sig2fdp50df), safe_mean(tempDat$sig2fdpNew),
          safe_mean(tempDat$sig2fdpDEIB), safe_mean(tempDat$sig2fdpMTEST)
        )
      }
    }
    
    out_list[[i]] <- summaryOut
  }
  
  dplyr::bind_rows(out_list)
}


# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# get output file names
outputDir <- here::here("Fig2", "output")

##output
names2a <- here::here(outputDir, paste0("Fig2A_aID", 120:780, ".txt"))
names2b <- here::here(outputDir, paste0("Fig2B_aID", 1:400, ".txt"))


#-----------------------------------------#

# read raw output files
res2a <- c()
for (file_it in 1:length(names2a)) {
  tempRes <- tryCatch(fread(names2a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res2a, tempRes)
  res2a <- rbindlist(tempList, fill = TRUE)
  # res2a <- na.omit(rbindlist(tempList, fill = TRUE))
}

res2b <- c()
for (file_it in 1:length(names2b)) {
  tempRes <- tryCatch(fread(names2b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res2b, tempRes)
  res2b <- rbindlist(tempList)
}


# summarize
summary2a <- summarize_raw(res2a)
summary2b <- summarize_raw(res2b)


# save summaries
write.table(summary2a, paste0(outputDir, "/Fig2a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary2b, paste0(outputDir, "/Fig2b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')


#-----------------------------------------#
# start plotting

# define colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mycols <- gg_color_hue(6)
mycols[4] <- "black"
mycols[5] <- "blue"

# Figure 2A
Fig2A_data <- fread(paste0(outputDir, "/Fig2a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  # filter(minEff1 >= 0.35 & minEff1 <= 0.59) %>%
  filter(!is.na(Method))

# plot Figure 2A
Fig2A_plot <- ggplot(data=Fig2A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (4D Mediation)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + #xlim(c(0.12, 0.22)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# plot Figure 2C
Fig2C_plot <- ggplot(data=Fig2A_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  # geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (4D Mediation)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 1.0)) + #xlim(c(0.12, 0.22)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))



# Figure 2B
Fig2B_data <- fread(paste0(outputDir, "/Fig2b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  # filter(minEff1 >= 0.35 & minEff1 <= 0.59) %>%
  filter(!is.na(Method))

# plot Figure 2B
Fig2B_plot <- ggplot(data=Fig2B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (4D Mediationn)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + #xlim(c(0.02, 0.11)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0,", delta[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# plot Figure 2D
Fig2D_plot <- ggplot(data=Fig2B_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  # geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (4D Mediation)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 1.0)) + #xlim(c(0.02, 0.11)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0,", delta[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))



# # put together figure 1
Fig2_plot <- plot_grid(Fig2A_plot + theme(legend.position = "none"),
                       Fig2B_plot + theme(legend.position = "none"),
                       Fig2C_plot + theme(legend.position = "none"),
                       Fig2D_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)
# 
# 
Fig2_legend <- get_legend(Fig2A_plot +  theme(legend.direction="horizontal",
                                              legend.justification="center",legend.box.just="bottom"))
# 
plot_grid(Fig2_plot, Fig2_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig2.pdf"), width=18, height=12)













