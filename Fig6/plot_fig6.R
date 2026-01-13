# Collect results and plot Figure 6
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/plot_fig3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again

setwd("~/Downloads/csmGmm_sim_R3/Fig6")
here::i_am("Fig6/plot_fig6.R")

library(ggplot2)
library(cowplot)
library(ggformula)
library(dplyr)
library(data.table)
library(csmGmm)
library(here)

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# get output file names
outputDir <- here::here("Fig6", "output")

##output
names6a <- here::here(outputDir, paste0("Fig6A_aID", 441:1041, ".txt"))
names6b <- here::here(outputDir, paste0("Fig6B_aID", 1:300, ".txt"))

#-----------------------------------------#

# read raw output files
res6a <- c()
for (file_it in 1:length(names6a)) {
  tempRes <- tryCatch(fread(names6a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res6a, tempRes)
  res6a <- rbindlist(tempList)
}

# Read 6b
res6b <- c()
for (file_it in 1:length(names6b)) {
  tempRes <- tryCatch(fread(names6b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res6b, tempRes)
  res6b <- rbindlist(tempList)
}


# summarize
summary6a <- summarize_raw(res6a)
summary6b <- summarize_raw(res6b)

# save summaries
write.table(summary6a, paste0(outputDir, "/Fig6a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary6b, paste0(outputDir, "/Fig6b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Figure 6A
Fig6A_data <- fread(paste0(outputDir, "/Fig6a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 >= 0.35 & minEff1 <= 0.59) %>%
  filter(!is.na(Method))

# plot Figure 6A
Fig6A_plot <- ggplot(data=Fig6A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Replication)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.3)) + xlim(c(0.36, 0.60)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# plot Figure 6C
Fig6C_plot <- ggplot(data=Fig6A_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Replication)") + 
  xlab("Min Effect Magnitude") +
  ylim(c(0, 1.0)) + xlim(c(0.36, 0.60)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

###

# Figure 6B or Correlation setting
Fig6B_data <- fread(paste0(outputDir, "/Fig6b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))


# plot Figure 6B
Fig6B_plot <- ggplot(data=Fig6B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Correlated Test Stat)") +
  ylim(c(0, 0.3)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# plot Figure 6D
Fig6D_plot <- ggplot(data=Fig6B_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Correlated Test Stat)") +
  ylim(c(0, 1.0)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# put together figure 1
Fig6_plot <- plot_grid(Fig6A_plot + theme(legend.position = "none"),
                       Fig6B_plot + theme(legend.position = "none"),
                       Fig6C_plot + theme(legend.position = "none"),
                       Fig6D_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)


Fig6_legend <- get_legend(Fig6A_plot +  theme(legend.direction="horizontal",
                                              legend.justification="center",legend.box.just="bottom"))

plot_grid(Fig6_plot, Fig6_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig6.pdf"), width=18, height=12)













