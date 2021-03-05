setwd("D:/SezermanLab/FLNA/cMDs/data/1_4m9p_cases/")

# RMSD ------------
rmsd_wt <- read.table("rmsd_4m9p.dat", sep = " ")
rmsd_l656f <- read.table("rmsd_L656F.dat", sep = " ")
rmsd_r484q <- read.table("rmsd_R484Q.dat", sep = " ")

rmsd_all <- cbind(rmsd_wt, rmsd_r484q[2], rmsd_l656f[2])
rmsd_all[1] <- seq(0, 99.99, 0.02)
colnames(rmsd_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")


rmsd_seperate <- reshape2::melt(rmsd_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rmsd_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 1.1, alpha = 0.75)
p <- p + geom_line(size = 1, alpha = 0.75)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 7, 0.5), limits = c(0,6.1))
p <- p + labs(x = "Time (ns)", y = "Deviations (Å)")
p <- p + ggsci::scale_color_startrek()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(size = 20),
               axis.title.x = element_text(size = 20),
               axis.ticks = element_line(),
               axis.text = element_text(size = 16),
               axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.position = c(0.9, 0.12),
               legend.key.width = unit(1.25,"cm"),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("RMSD_4m9pCases.pdf", width = 10, height = 6)
p
dev.off()

# RMSF -------------
rmsf_wt <- read.table("rmsf_4m9p.dat", sep = " ")
rmsf_l656f <- read.table("rmsf_L656F.dat", sep = " ")
rmsf_r484q <- read.table("rmsf_R484Q.dat", sep = " ")

rmsf_all <- cbind(rmsf_wt, rmsf_r484q[2], rmsf_l656f[2])
rmsf_all[1] <- c(478:766)
colnames(rmsf_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

rmsf_seperate <- reshape2::melt(rmsf_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rmsf_seperate, aes(x=time, y=value, colour=variable))
p <- p + geom_point(size = 1.1, alpha = 0.75)
p <- p + geom_line(size = 1, alpha = 0.75)
p <- p + scale_x_continuous(breaks = c(478,seq(490,750,20),766), limits = c(472, 772),
                            expand = c(0,0), minor_breaks = 472:772)
p <- p + scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 9), expand = c(0,0.1))
p <- p + labs(x = "Residue Index", y = "Fluctuations (Å)")
p <- p + ggsci::scale_color_startrek()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 20),
               axis.title.x = element_text(size = 20),
               axis.ticks = element_line(),
               axis.text = element_text(size = 16),
               axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.position = c(0.88, 0.92),
               legend.key.width = unit(1.25,"cm"),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("RMSF_4m9pCases.pdf", width = 10, height = 6)
p
dev.off()

# RMSF (focused) -------------
rmsf_wt <- read.table("rmsf_4m9p.dat", sep = " ")
rmsf_l656f <- read.table("rmsf_L656F.dat", sep = " ")
rmsf_r484q <- read.table("rmsf_R484Q.dat", sep = " ")

rmsf_all <- cbind(rmsf_wt, rmsf_r484q[2], rmsf_l656f[2])
rmsf_all[1] <- c(478:766)
colnames(rmsf_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

rmsf_seperate <- reshape2::melt(rmsf_all[1:4], id.var = "time")
library(ggplot2)
library(ggpubr)
p1 <- ggplot(rmsf_seperate, aes(x=time, y=value, colour=variable))
p1 <- p1 + geom_point(size = 1.1, alpha = 0.75)
p1 <- p1 + geom_line(size = 1, alpha = 0.75)
p1 <- p1 + scale_x_continuous(breaks = seq(505,570,5), limits = c(500,575),
                              expand = c(0,0), minor_breaks = 500:575)
p1 <- p1 + scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 5), expand = c(0,0.1))
p1 <- p1 + labs(x = "Residue Index", y = "Fluctuations (Å)")
p1 <- p1 + ggsci::scale_color_startrek()
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(#axis.title.y = element_text(angle = 90, size = 20),
                 # axis.title.x = element_text(size = 20),
                 axis.title = element_blank(),
                 axis.ticks = element_line(),
                 axis.text = element_text(size = 16),
                 axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 # legend.position = c(0.88, 0.92),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1

p2 <- ggplot(rmsf_seperate, aes(x=time, y=value, colour=variable))
p2 <- p2 + geom_point(size = 1.1, alpha = 0.75)
p2 <- p2 + geom_line(size = 1, alpha = 0.75)
p2 <- p2 + scale_x_continuous(breaks = c(624,630), limits = c(623,631),
                              expand = c(0,0), minor_breaks = 623:631)
p2 <- p2 + scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 5), expand = c(0,0.1))
p2 <- p2 + labs(x = "Residue Index", y = "Fluctuations (Å)")
p2 <- p2 + ggsci::scale_color_startrek()
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(#axis.title.y = element_text(angle = 90, size = 20),
                 # axis.title.x = element_text(size = 20),
                 axis.title = element_blank(),
                 axis.ticks.x = element_line(),
                 axis.ticks.y = element_line(),
                 axis.text.y = element_blank(),
                 axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 16),
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12),
                 # legend.position = c(0.88, 0.92),
                 legend.position = "none",
                 legend.key.width = unit(1.25,"cm"),
                 panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2

p3 <- ggplot(rmsf_seperate, aes(x=time, y=value, colour=variable))
p3 <- p3 + geom_point(size = 1.1, alpha = 0.75)
p3 <- p3 + geom_line(size = 1, alpha = 0.75)
p3 <- p3 + scale_x_continuous(breaks = seq(685,740,5), limits = c(685,742),
                              expand = c(0,0), minor_breaks = 685:742)
p3 <- p3 + scale_y_continuous(breaks = seq(0, 10, 1), limits = c(0, 5), expand = c(0,0.1))
p3 <- p3 + labs(x = "Residue Index", y = "Fluctuations (Å)")
p3 <- p3 + ggsci::scale_color_startrek()
p3 <- p3 + theme_minimal()
p3 <- p3 + theme(#axis.title.y = element_text(angle = 90, size = 20),
  # axis.title.x = element_text(size = 20),
  axis.title = element_blank(),
  axis.ticks.x = element_line(),
  axis.ticks.y = element_line(),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 16),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = c(0.77, 0.92),
  legend.key.width = unit(1.25,"cm"),
  panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p3

fig <- ggarrange(p1,p2,p3,widths = c(5.357,0.571,4.071), nrow = 1, ncol = 3, align = "h")

pdf("RMSF_4m9pCases_focused_part1-3.pdf", width = 10, height = 6)
annotate_figure(fig, left = text_grob("Fluctuations (Å)", size = 20, rot = 90),
                bottom = text_grob("Residue Index", size = 20))
dev.off()

# ROG --------------
rog_wt <- read.table("rog_4m9p.dat", sep = "\t")
rog_l656f <- read.table("rog_L656F.dat", sep = "\t")
rog_r484q <- read.table("rog_R484Q.dat", sep = "\t")

rog_all <- cbind(rog_wt, rog_r484q[2], rog_l656f[2])
rog_all[1] <- seq(0, 99.99, 0.02)
colnames(rog_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

library(reshape2)
rog_seperate <- melt(rog_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rog_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 1.1, alpha = 0.75)
p <- p + geom_line(size = 1, alpha = 0.75)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 50, 0.5), limits = c(19.5,22.25))
p <- p + labs(x = "Time (ns)", y = "Compactness (Å)")
p <- p + ggsci::scale_color_startrek()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(size = 20),
               axis.title.x = element_text(size = 20),
               axis.ticks = element_line(),
               axis.text = element_text(size = 16),
               axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.position = c(0.1, 0.12),
               legend.key.width = unit(1.25,"cm"),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("Rg_4m9pCases.pdf", width = 10, height = 6)
p
dev.off()

# stability --------------
stability_wt <- read.table("stability_4m9p.dat", sep = "\t")
stability_l656f <- read.table("stability_L656F.dat", sep = "\t")
stability_r484q <- read.table("stability_R484Q.dat", sep = "\t")

stability_all <- cbind(stability_wt, stability_r484q[2], stability_l656f[2])
stability_all[1] <- seq(0, 99.99, 1.960784)
colnames(stability_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

library(reshape2)
stability_seperate <- melt(stability_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(stability_seperate, aes(x=variable, y=value, colour=variable))
p <- p + geom_boxplot(size = 1, alpha = 0.75)
p <- p + scale_y_continuous(breaks = seq(0, 1500, 5), limits = c(110,165))
p <- p + labs(x = "Time (ns)", y = expression("FoldX Stability ("*Delta*"G, kcal/mol)"))
p <- p + ggsci::scale_color_startrek()
# p <- p + stat_compare_means(comparisons = comp, method = "wilcox.test", step.increase = 0.14)
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(size = 20),
               axis.title.x = element_blank(),
               axis.ticks = element_line(),
               axis.text = element_text(size = 16),
               axis.text.x = element_text(size = 13, hjust = 0.5, vjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.position = "none",
               legend.key.width = unit(1.25,"cm"),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("stability_4m9pCases.pdf", width = 5, height = 6)
p
dev.off()

# SASA --------------
SASA_wt <- read.table("SASA_overall_4m9p.dat", sep = "\t")
SASA_l656f <- read.table("SASA_overall_L656F.dat", sep = "\t")
SASA_r484q <- read.table("SASA_overall_R484Q.dat", sep = "\t")

SASA_all <- cbind("",SASA_wt, SASA_r484q, SASA_l656f)
SASA_all[1] <- seq(0, 99.99, 0.02)
colnames(SASA_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

library(reshape2)
SASA_seperate <- melt(SASA_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(SASA_seperate, aes(x=variable, y=value, colour=variable))
p <- p + geom_boxplot(size = 1, alpha = 0.75)
p <- p + scale_y_continuous(breaks = seq(14500, 16750, 250), limits = c(14500,16750))
p <- p + labs(x = "Time (ns)", y = "SASA (Å²)")
p <- p + ggsci::scale_color_startrek()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(size = 20),
               axis.title.x = element_blank(),
               axis.ticks = element_line(),
               axis.text = element_text(size = 16),
               axis.text.x = element_text(size = 13, hjust = 0.5, vjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 12),
               legend.position = "none",
               legend.key.width = unit(1.25,"cm"),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("SASA_4m9pCases.pdf", width = 5, height = 6)
p
dev.off()