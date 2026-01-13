# Collect results and plot Figure 5
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
Fig5A_data <- fread(paste0(outputDir, "/Fig3a_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.25 & minEff1 <= 0.475) %>%
  filter(!is.na(Method))

# plot Figure 
Fig5A_plot_incon <- ggplot(data=Fig5A_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Num Incongruous (3D Pleiotropy)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 3900)) + xlim(c(0.26, 0.48)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# Figure 5B
Fig3B_data <- fread(paste0(outputDir, "/Fig3b_summary.txt"), data.table=F) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))

Fig5B_plot_incon <- ggplot(data=Fig5B_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Num Incongruous (3D Pleiotropy)") +
  ylim(c(0, 11500)) +
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))


# Figure 5C
Fig5C_data <-  read.table("~/Downloads/csmGmm_sim_R3/Fig6/output/Fig6a_summary.txt", header = TRUE) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  filter(minEff1 >= 0.25 & minEff1 <= 0.475) %>%
  filter(!is.na(Method))


# plot Figure 
Fig5C_plot_incon <- ggplot(data=Fig5C_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) +
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Num Incongruous (3D Replication)") +
  xlab("Min Effect Magnitude") +
  ylim(c(0, 2000)) + xlim(c(0.305, 0.48)) +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18)) +
  theme(legend.key.size = unit(3,"line"))


# Figure 5D
Fig5D_data <- read.table("~/Downloads/csmGmm_sim_R3/Fig6/output/Fig6b_summary.txt", header = TRUE) %>%
  mutate(Method = ifelse(Method == "df7", "locfdr7df", Method)) %>%
  mutate(Method = ifelse(Method == "df50", "locfdr50df", Method)) %>%
  mutate(Method = ifelse(Method == "New", "csmGmm", Method)) %>%
  mutate(Method = factor(Method, levels=c("csmGmm", "Kernel", "locfdr7df",
                                          "locfdr50df")))  %>%
  #filter(minEff1 <= 0.08) %>%
  filter(!is.na(Method))

Fig5D_plot_incon <- ggplot(data=Fig5D_data, aes(x=minEff1, y=Incongruous, group=Method)) +
  geom_line(aes(linetype = Method, color=Method),lwd=1.2) + 
  scale_color_manual(values=mycols[-c(2,6)]) +
  scale_linetype_manual(values=c(1,2,3,4,5,6)[-c(2,6)]) +
  #geom_hline(yintercept=0.1, linetype=2, color="grey") +
  ylab("Num Incongruous (3D Correlated Test Stat)") +
  ylim(c(0, 4000)) +
  xlab("Min Effect Magnitude") +
  theme_cowplot() +
  theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=18))+
  theme(legend.key.size = unit(3,"line"))


# put together figure 5
Fig5_plot <- plot_grid(Fig5A_plot_incon + theme(legend.position = "none"),
                       Fig5B_plot_incon + theme(legend.position = "none"),
                       Fig5C_plot_incon + theme(legend.position = "none"),
                       Fig5D_plot_incon + theme(legend.position = "none"),
                       labels=c("A", "B", "C", "D"), nrow=2, label_size=22)


Fig5_legend <- get_legend(Fig5A_plot_incon + theme(legend.direction="horizontal",
                                             legend.justification="center",legend.box.just="bottom"))

plot_grid(Fig5_plot, Fig5_legend, ncol=1, rel_heights=c(1, 0.1))
ggsave(paste0(outputDir, "/Fig5.pdf"), width=18, height=12)



