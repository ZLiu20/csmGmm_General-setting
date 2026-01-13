# Collect results and plot Figure 3
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/plot_fig3.R") or set the path after the -cwd flag
# in the .lsf file, and then run again.

setwd("~/Downloads/csmGmm_sim_R3/Fig3")
here::i_am("Fig3/plot_fig3.R")

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
outputDir <- here::here("Fig3", "output")

##output
names1a <- here::here(outputDir, paste0("Fig3A_aID", 221:720, ".txt"))
names1b <- here::here(outputDir, paste0("Fig3B_aID", 1:200, ".txt"))

#-----------------------------------------#

# read raw output files
res3a <- c()
for (file_it in 1:length(names3a)) {
  tempRes <- tryCatch(fread(names3a[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3a, tempRes)
  res3a <- rbindlist(tempList)
}

# Read 3b
res3b <- c()
for (file_it in 1:length(names3b)) {
  tempRes <- tryCatch(fread(names3b[file_it]), error=function(e) e)
  if (class(tempRes)[1] == "simpleError") {next}
  tempList <- list(res3b, tempRes)
  res3b <- rbindlist(tempList)
}


# summarize
summary3a <- summarize_raw(res3a)
summary3b <- summarize_raw(res3b)

# save summaries
write.table(summary3a, paste0(outputDir, "/Fig3a_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')
write.table(summary3b, paste0(outputDir, "/Fig3b_summary.txt"), append=F,quote=F, row.names=F, col.names=T, sep='\t')

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

# Figure 3A
Fig3A_data <- fread(paste0(outputDir, "/Fig3a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.25 & minEff1 <= 0.475) %>%
  filter(!is.na(Method))

# plot Figure 3A
Fig3A_plot <- ggplot(data=Fig3A_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 0.2)) + xlim(c(0.26, 0.48)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# plot Figure 3C
Fig3C_plot <- ggplot(data=Fig3A_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 1.0)) + xlim(c(0.26, 0.48)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))

# Figure 3B
Fig3B_data <- fread(paste0(outputDir, "/Fig3b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))

# plot Figure 3B
Fig3B_plot <- ggplot(data=Fig3B_data, aes(x=minEff1, y=FDP, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("FDP (3D Pleiotropy)") +
  ylim(c(0, 0.2)) + #xlim(c(0.01, 0.03 )) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# plot Figure 3D
Fig3D_plot <- ggplot(data=Fig3B_data, aes(x=minEff1, y=Power, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Power (3D Pleiotropy)") +
  ylim(c(0, 1.0)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))

# put together figure 1
Fig1_plot <- plot_grid(Fig3A_plot + theme(legend.position = "none"),
                        Fig3B_plot + theme(legend.position = "none"),
                        Fig3C_plot + theme(legend.position = "none"),
                        Fig3D_plot + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)


Fig3_legend <- get_legend(Fig3A_plot +  theme(legend.direction="horizontal",
                                              legend.justification="center",legend.box.just="bottom"))

plot_grid(Fig3_plot, Fig3_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig3.pdf"), width=18, height=12)













