# Collect results and plot Figure 4
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/plot_fig3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.

setwd("~/Downloads/csmGmm_sim_R3/Fig4")
here::i_am("Fig4/plot_fig4.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(csmGmm)
library(here)



# Function to summarize raw results
summarize_raw_R1 <- function(fullDat, full=FALSE, cor=FALSE, maxP=FALSE, FDP2=FALSE) {
  # summarize raw output results
  outDF <- c()
  effSizes <- sort(unique(fullDat$minEff1))
  for (eff_it in 1:length(effSizes)) {
    # loop through each effect size
    tempEff <- effSizes[eff_it]
    tempDat <- fullDat %>% filter(minEff1 == tempEff)  %>%
      as.data.frame(.) %>%
      # mutate(fdpDACT = ifelse(nRejDACT == 0, 0, fdpDACT)) %>%
      # mutate(fdpHDMT = ifelse(nRejHDMT == 0, 0, fdpHDMT)) %>%
      mutate(fdpKernel = ifelse(nRejKernel == 0, 0, fdpKernel)) %>%
      mutate(fdp7df = ifelse(nRej7df == 0, 0, fdp7df)) %>%
      mutate(fdp50df = ifelse(nRej50df == 0, 0, fdp50df)) %>%
      mutate(fdpNew = ifelse(nRejNew == 0, 0, fdpNew))
    
    if (cor) {
      tempDat <- tempDat %>% mutate(fdpCor = ifelse(nRejCor == 0, 0, fdpCor)) %>%
        mutate(fdpNew = fdpCor) %>%
        mutate(powerNew = powerCor) %>%
        mutate(nRejNew = nRejCor)
    }
    
    if (full) {
      tempDat <- tempDat %>% mutate(fdpFull = ifelse(nRejFull == 0, 0, fdpFull)) %>%
        mutate(fdpDACT = fdpFull) %>%
        mutate(powerDACT = powerFull) %>%
        mutate(nRejDACT = nRejFull)
    }
    
    # summarize
    summaryOut <- data.frame(minEff1 = tempDat$minEff1[1],
                             #Method = c("DACT", "HDMT", "Kernel", "df7", "df50", "New"))
                             Method = c("Kernel", "df7","df50","New"))
    summaryOut$nRej <- c(#mean(tempDat$nRejDACT, na.rm=T), mean(tempDat$nRejHDMT, na.rm=T),
                         mean(tempDat$nRejKernel, na.rm=T), mean(tempDat$nRej7df, na.rm=T),
                         mean(tempDat$nRej50df, na.rm=T),
                         mean(tempDat$nRejNew, na.rm=T))
    summaryOut$Power <- c(#mean(tempDat$powerDACT, na.rm=T), mean(tempDat$powerHDMT, na.rm=T),
                        mean(tempDat$powerKernel, na.rm=T), mean(tempDat$power7df, na.rm=T), 
                        mean(tempDat$power50df, na.rm=T),
                          mean(tempDat$powerNew, na.rm=T))
    summaryOut$FDP <- c(#mean(tempDat$fdpDACT, na.rm=T),
                        # mean(tempDat$fdpHDMT, na.rm=T), 
                        mean(tempDat$fdpKernel, na.rm=T), mean(tempDat$fdp7df, na.rm=T), 
                        mean(tempDat$fdp50df, na.rm=T),
                        mean(tempDat$fdpNew, na.rm=T))
    summaryOut$Incongruous <- c(#NA, NA, 
                                mean(tempDat$inconKernel, na.rm=T), mean(tempDat$incon7df, na.rm=T),
                                mean(tempDat$incon50df, na.rm=T), 
                                mean(tempDat$inconNew, na.rm=T))
    summaryOut$numNA <- c(#length(which(is.na(tempDat$powerDACT))),
                          # length(which(is.na(tempDat$powerHDMT))), 
                          length(which(is.na(tempDat$powerKernel))),
                          length(which(is.na(tempDat$power7df))), 
                          length(which(is.na(tempDat$power50df))), 
                          length(which(is.na(tempDat$powerNew))))
    
    if (maxP) {
      summaryOut <- data.frame(minEff1 = tempDat$minEff1[1],
                               #Method = c("MaxP", "DACT", "HDMT", "Kernel", "df7", "df50", "New"))
                               Method = c("MaxP", "DACT", "HDMT", "Kernel", "df7", "df50", "New"))
      summaryOut$nRej <- c(mean(tempDat$nRejMaxp, na.rm=T), mean(tempDat$nRejDACT, na.rm=T), mean(tempDat$nRejHDMT, na.rm=T),
                           mean(tempDat$nRejKernel, na.rm=T), mean(tempDat$nRej7df, na.rm=T),
                           mean(tempDat$nRej50df, na.rm=T), mean(tempDat$nRejNew, na.rm=T))
      summaryOut$Power <- c(mean(tempDat$powerMaxp, na.rm=T), mean(tempDat$powerDACT, na.rm=T),
                            mean(tempDat$powerHDMT, na.rm=T), mean(tempDat$powerKernel, na.rm=T),
                            mean(tempDat$power7df, na.rm=T), mean(tempDat$power50df, na.rm=T),
                            mean(tempDat$powerNew, na.rm=T))
      summaryOut$FDP <- c(mean(tempDat$fdpMaxp, na.rm=T), mean(tempDat$fdpDACT, na.rm=T),
                          mean(tempDat$fdpHDMT, na.rm=T), mean(tempDat$fdpKernel, na.rm=T),
                          mean(tempDat$fdp7df, na.rm=T), mean(tempDat$fdp50df, na.rm=T),
                          mean(tempDat$fdpNew, na.rm=T))
    }
    
    if (FDP2) {
      summaryOut$sig2Pow <- c(mean(tempDat$sig2PowDACT, na.rm=T),
                              mean(tempDat$sig2PowHDMT, na.rm=T), mean(tempDat$sig2PowKernel, na.rm=T),
                              mean(tempDat$sig2Pow7df, na.rm=T), mean(tempDat$sig2Pow50df, na.rm=T),
                              mean(tempDat$sig2PowNew, na.rm=T))
      summaryOut$sig2FDP <- c(mean(tempDat$sig2fdpDACT, na.rm=T),
                              mean(tempDat$sig2fdpHDMT, na.rm=T), mean(tempDat$sig2fdpKernel, na.rm=T),
                              mean(tempDat$sig2fdp7df, na.rm=T), mean(tempDat$sig2fdp50df, na.rm=T),
                              mean(tempDat$sig2fdpNew, na.rm=T))
    }
    outDF <- rbind(outDF, summaryOut)
  }
  
  return(outDF)
}



# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# get output file names
outputDir <- here::here("Fig4", "output")

##output
names4a <- here::here(outputDir, paste0("Fig4A_aID", 1:1000, ".txt"))
names4b <- here::here(outputDir, paste0("Fig4B_aID", 1:1000, ".txt"))

#-----------------------------------------#

# read raw output files

res4a <- list() 
for (file_it in seq_along(names4a)) {
  tempRes <- tryCatch(fread(names4a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  
  if (is.data.frame(tempRes)) {  
    res1a[[file_it]] <- tempRes  
  }
}

res4a <- rbindlist(res4a, use.names=TRUE, fill=TRUE)  


res4b <- list() 
for (file_it in seq_along(names4b)) {
  tempRes <- tryCatch(fread(names1b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  
  if (is.data.frame(tempRes)) {  
    res4b[[file_it]] <- tempRes  
  }
}

res4b <- rbindlist(res4b, use.names=TRUE, fill=TRUE) 


res4a <- res4a %>% filter(!is.na(res4a$power50df) & !is.na(res1a$fdp50df) )
res4b <- res4b %>% filter(!is.na(res4b$power50df) & !is.na(res1b$fdp50df) )


# summarize
summary4a <- summarize_raw_R1(res4a)
summary4b <- summarize_raw_R1(res4b)


# save summaries
write.table(summary4a, paste0(outputDir, "/Fig4a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary4b, paste0(outputDir, "/Fig4b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Figure 4A
Fig4A_data <- fread(paste0(outputDir, "/Fig4a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.35 & minEff1 <= 0.75) %>%
  filter(!is.na(Method))

# plot Figure 4A
Fig4A_plot <- ggplot(data=Fig4A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (4D Pleiotropy)") +
  xlab("Min Effect Magnitude (t=1)") +
  ylim(c(0, 0.3)) + xlim(c(0.35, 0.70)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# plot Figure 4C
Fig4C_plot <- ggplot(data=Fig5A_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (4D Pleiotropy)") +
  xlab("Min Effect Magnitude (t=1)") +
  ylim(c(0, 1.0)) + xlim(c(0.35, 0.70)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# Figure 4B
Fig4B_data <- fread(paste0(outputDir, "/Fig4b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.35 & minEff1 <= 0.75) %>%
  filter(!is.na(Method))

# plot Figure 4B
Fig4B_plot <- ggplot(data=Fig4B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (4D Pleiotropy)") +
  xlab("Min Effect Magnitude (t=2)") +
  ylim(c(0, 0.45)) + xlim(c(0.35, 0.70)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# plot Figure 4D
Fig4D_plot <- ggplot(data=Fig4B_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  # geom_hline(yintercept=0.1, linetype=2, color="grey") +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (4D Pleiotropy)") +
  xlab("Min Effect Magnitude (t=2)") +
  ylim(c(0, 1.0)) + xlim(c(0.35, 0.70)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))


# put together figure 1
Fig4_plot <- plot_grid(Fig4A_plot + theme(legend.position = "none"),
                       Fig4B_plot + theme(legend.position = "none"),
                       Fig4C_plot + theme(legend.position = "none"),
                       Fig4D_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)


Fig4_legend <- get_legend(Fig4A_plot +  theme(legend.direction="horizontal",
                                              legend.justification="center",legend.box.just="bottom"))

plot_grid(Fig4_plot, Fig4_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig4.pdf"), width=18, height=12)













