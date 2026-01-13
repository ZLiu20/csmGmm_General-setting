# Collect results and plot Figure 3
# Using the here package to manage file paths. If an error is thrown, please
# set the working directory to the folder that holds this Rscript, e.g.
# setwd("/path/to/csmGmm_reproduce/Fig3/plot_fig3.R") or set the path after the -cwd flag
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

# source the .R scripts from the SupportingCode/ folder 
codePath <- c(here::here("SupportingCode"))
toBeSourced <- list.files(codePath, "\\.R$")
purrr::map(paste0(codePath, "/", toBeSourced), source)

# get output file names
outputDir <- here::here("Fig2", "output")

##output
names2a <- here::here(outputDir, paste0("Fig7A_aID", 1:841, ".txt"))
names2b <- here::here(outputDir, paste0("Fig7B_aID", 1:400, ".txt"))


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
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
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
  xlab(expression(paste(tau[1] ,"= Proportion of (", alpha[j] != 0, ", ", beta[j], "=0, ", gamma[j], "=0)"))) +
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













