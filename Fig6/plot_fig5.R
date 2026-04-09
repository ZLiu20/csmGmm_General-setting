# Collect results and plot Figure 5
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/plot_fig3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again

setwd("~/Downloads/csmGmm_sim_R3/Fig5")
here::i_am("Fig5/plot_fig5.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(csmGmm)
library(here)



# Function to summarize raw results
summarize_raw_R1 <- function(fullDat, full = FALSE, cor = FALSE, maxP = FALSE, FDP2 = FALSE) {
  outDF <- data.frame()
  effSizes <- sort(unique(fullDat$minEff1))
  
  needed <- c(
    "nRejDACT","nRejHDMT","nRejKernel","nRej7df","nRej50df","nRejNew","nRejQCH","nRejAda",
    "powerDACT","powerHDMT","powerKernel","power7df","power50df","powerNew","powerQCH","powerAda",
    "fdpDACT","fdpHDMT","fdpKernel","fdp7df","fdp50df","fdpNew","fdpQCH","fdpAda",
    "inconKernel","incon7df","incon50df","inconNew"
  )
  for (nm in needed) if (!nm %in% names(fullDat)) fullDat[[nm]] <- NA_real_
  
  for (tempEff in effSizes) {
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
        fdpQCH    = ifelse(nRejQCH    == 0, 0, fdpQCH),
        fdpAda    = ifelse(nRejAda    == 0, 0, fdpAda)
      )
    
    methods <- if (maxP) c("MaxP","DACT","HDMT","Kernel","df7","df50","New","QCH","Ada") else
      c("DACT","HDMT","Kernel","df7","df50","New","QCH","Ada")
    
    summaryOut <- data.frame(minEff1 = tempDat$minEff1[1], Method = methods)
    
    if (maxP) {
      summaryOut$nRej <- c(mean(tempDat$nRejMaxp, na.rm=TRUE), mean(tempDat$nRejDACT, na.rm=TRUE), mean(tempDat$nRejHDMT, na.rm=TRUE),
                           mean(tempDat$nRejKernel, na.rm=TRUE), mean(tempDat$nRej7df, na.rm=TRUE), mean(tempDat$nRej50df, na.rm=TRUE),
                           mean(tempDat$nRejNew, na.rm=TRUE), mean(tempDat$nRejQCH, na.rm=TRUE), mean(tempDat$nRejAda, na.rm=TRUE))
      summaryOut$Power <- c(mean(tempDat$powerMaxp, na.rm=TRUE), mean(tempDat$powerDACT, na.rm=TRUE), mean(tempDat$powerHDMT, na.rm=TRUE),
                            mean(tempDat$powerKernel, na.rm=TRUE), mean(tempDat$power7df, na.rm=TRUE), mean(tempDat$power50df, na.rm=TRUE),
                            mean(tempDat$powerNew, na.rm=TRUE), mean(tempDat$powerQCH, na.rm=TRUE), mean(tempDat$powerAda, na.rm=TRUE))
      summaryOut$FDP <- c(mean(tempDat$fdpMaxp, na.rm=TRUE), mean(tempDat$fdpDACT, na.rm=TRUE), mean(tempDat$fdpHDMT, na.rm=TRUE),
                          mean(tempDat$fdpKernel, na.rm=TRUE), mean(tempDat$fdp7df, na.rm=TRUE), mean(tempDat$fdp50df, na.rm=TRUE),
                          mean(tempDat$fdpNew, na.rm=TRUE), mean(tempDat$fdpQCH, na.rm=TRUE), mean(tempDat$fdpAda, na.rm=TRUE))
      summaryOut$Incongruous <- c(NA, NA, NA, mean(tempDat$inconKernel, na.rm=TRUE), mean(tempDat$incon7df, na.rm=TRUE),
                                  mean(tempDat$incon50df, na.rm=TRUE), mean(tempDat$inconNew, na.rm=TRUE), 
                                  mean(tempDat$inconQCH, na.rm=TRUE), mean(tempDat$inconAda, na.rm=TRUE))
      summaryOut$numNA <- c(length(which(is.na(tempDat$powerMaxp))), length(which(is.na(tempDat$powerDACT))), length(which(is.na(tempDat$powerHDMT))),
                            length(which(is.na(tempDat$powerKernel))), length(which(is.na(tempDat$power7df))), length(which(is.na(tempDat$power50df))),
                            length(which(is.na(tempDat$powerNew))), length(which(is.na(tempDat$powerQCH))), length(which(is.na(tempDat$powerAda))))
    } else {
      summaryOut$nRej <- c(mean(tempDat$nRejDACT, na.rm=TRUE), mean(tempDat$nRejHDMT, na.rm=TRUE), mean(tempDat$nRejKernel, na.rm=TRUE),
                           mean(tempDat$nRej7df, na.rm=TRUE), mean(tempDat$nRej50df, na.rm=TRUE), mean(tempDat$nRejNew, na.rm=TRUE),
                           mean(tempDat$nRejQCH, na.rm=TRUE), mean(tempDat$nRejAda, na.rm=TRUE))
      summaryOut$Power <- c(mean(tempDat$powerDACT, na.rm=TRUE), mean(tempDat$powerHDMT, na.rm=TRUE), mean(tempDat$powerKernel, na.rm=TRUE),
                            mean(tempDat$power7df, na.rm=TRUE), mean(tempDat$power50df, na.rm=TRUE), mean(tempDat$powerNew, na.rm=TRUE),
                            mean(tempDat$powerQCH, na.rm=TRUE), mean(tempDat$powerAda, na.rm=TRUE))
      summaryOut$FDP <- c(mean(tempDat$fdpDACT, na.rm=TRUE), mean(tempDat$fdpHDMT, na.rm=TRUE), mean(tempDat$fdpKernel, na.rm=TRUE),
                          mean(tempDat$fdp7df, na.rm=TRUE), mean(tempDat$fdp50df, na.rm=TRUE), mean(tempDat$fdpNew, na.rm=TRUE),
                          mean(tempDat$fdpQCH, na.rm=TRUE), mean(tempDat$fdpAda, na.rm=TRUE))
      summaryOut$Incongruous <- c(NA, NA, mean(tempDat$inconKernel, na.rm=TRUE), mean(tempDat$incon7df, na.rm=TRUE),
                                  mean(tempDat$incon50df, na.rm=TRUE), mean(tempDat$inconNew, na.rm=TRUE), 
                                  mean(tempDat$inconQCH, na.rm=TRUE), mean(tempDat$inconAda, na.rm=TRUE))
      summaryOut$numNA <- c(length(which(is.na(tempDat$powerDACT))), length(which(is.na(tempDat$powerHDMT))), length(which(is.na(tempDat$powerKernel))),
                            length(which(is.na(tempDat$power7df))), length(which(is.na(tempDat$power50df))), length(which(is.na(tempDat$powerNew))),
                            length(which(is.na(tempDat$powerQCH))), length(which(is.na(tempDat$powerAda))))
    }
    
    outDF <- rbind(outDF, summaryOut)
  }
  
  outDF
}



# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# get output file names
outputDir <- here::here("Fig5", "output")

##output
names5a <- here::here(outputDir, paste0("Fig5a_aID", 441:1041, ".txt"))
names5b <- here::here(outputDir, paste0("Fig5b_aID", 310:1041, ".txt"))

#-----------------------------------------#

# read raw output files
res5a <- c()
for (file_it in 1:length(names5a)) {
  tempRes <- tryCatch(fread(names5a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res5a, tempRes)
  res5a <- rbindlist(tempList)
}

# Read 5b
res5b <- c()
for (file_it in 1:length(names5b)) {
  tempRes <- tryCatch(fread(names5b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res5b, tempRes)
  res5b <- rbindlist(tempList)
}


# summarize
summary5a <- summarize_raw(res5a)
summary5b <- summarize_raw_R1(res5b)

# save summaries
write.table(summary5a, paste0(outputDir, "/Fig5a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary5b, paste0(outputDir, "/Fig5b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Figure 5A
Fig5A_data <- fread(paste0(outputDir, "/Fig5a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 >= 0.35 & minEff1 <= 0.59) %>%
  filter(!is.na(Method))

# plot Figure 5A
Fig5A_plot <- ggplot(data=Fig5A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Replication)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + #xlim(c(0.36, 0.60)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# plot Figure 5C
Fig5C_plot <- ggplot(data=Fig5A_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Replication)") + 
  xlab("Min Effect Magnitude") +
  ylim(c(0, 1.0)) + #xlim(c(0.36, 0.60)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# plot Figure  5E
Fig5E_plot <- ggplot(data=Fig5A_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Num Incongruous (3D Replication)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 2000)) + #xlim(c(0.305, 0.48)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

###

# Figure 5B or Correlation setting
Fig5B_data <- fread(paste0(outputDir, "/Fig5b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = ifelse(Method == "QCH", "qch_copula", Method)) %>%
  mutate(Method = ifelse(Method == "Ada", "adaFilter", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "adaFilter", "Kernel", "locfdr7df", "locfdr50df", "qch_copula"))) %>%
  #filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))


# plot Figure 5B
Fig5B_plot <- ggplot(data=Fig5B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Correlated Test Stat)") +
  ylim(c(0, 0.35)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# plot Figure 5D
Fig5D_plot <- ggplot(data=Fig5B_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Correlated Test Stat)") +
  ylim(c(0, 1.0)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# plot Figure 5F
Fig5F_plot <- ggplot(data=Fig5B_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  # scale_color_manual(values=mycols[-c(2,6)]) +
  # scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  scale_color_manual(values=mycols) +
  scale_linetype_manual(values=1:6) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Correlated Test Stat)") +
  ylim(c(0, 4000)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

                      
# put together figure 1
Fig5_plot <- plot_grid(Fig5A_plot + theme(legend.position = "none"),
                       Fig5B_plot + theme(legend.position = "none"),
                       Fig5C_plot + theme(legend.position = "none"),
                       Fig5D_plot + theme(legend.position = "none"),
                       Fig5E_plot + theme(legend.position = "none"),
                       Fig5F_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D", "E", "F"), nrow=3, label_size=22)


Fig5_legend <- get_legend(Fig5A_plot +  theme(legend.direction="horizontal",
                                              legend.justification="center",legend.box.just="bottom"))

plot_grid(Fig5_plot, Fig5_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig5.pdf"), width=14, height=18)













