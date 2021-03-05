setwd("D:/SezermanLab/FLNA/cMDs/data/1_4m9p_cases/")

# SASA (overall) ------------------
SASA_wt <- read.table("SASA_overall_4m9p.dat", sep = " ")
SASA_l656f <- read.table("SASA_overall_L656F.dat", sep = " ")
SASA_r484q <- read.table("SASA_overall_R484Q.dat", sep = " ")

SASA_all <- cbind("", SASA_wt, SASA_l656f, SASA_r484q)
SASA_all[1] <- seq(0, 99.99, 0.02)
colnames(SASA_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

SASA_seperate <- reshape2::melt(SASA_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(SASA_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.75) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 7, 0.5))
p <- p + labs(x = "Time (ns)", y = "Solvent Accessible Surface Area (Å²)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 21),
               axis.title.x = element_text(size = 21),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, 1.01),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("SASA_overall_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# SASA btw IgFLNa3-4 -----------------
SASA_wt <- read.table("SASA_btwIgFLNa3_4_4m9p.dat", sep = " ")
SASA_l656f <- read.table("SASA_btwIgFLNa3_4_L656F.dat", sep = " ")
SASA_r484q <- read.table("SASA_btwIgFLNa3_4_R484Q.dat", sep = " ")

SASA_all <- cbind("", SASA_wt, SASA_l656f, SASA_r484q)
SASA_all[1] <- seq(0, 99.99, 0.02)
colnames(SASA_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

SASA_seperate <- reshape2::melt(SASA_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(SASA_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.75) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 7, 0.5))
p <- p + labs(x = "Time (ns)", y = "Solvent Accessible Surface Area (Å²)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 21),
               axis.title.x = element_text(size = 21),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, 1.01),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("SASA_btwIgFLNa3-4_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Rg IgFLNa3 ---------------
rog_wt <- read.table("Rg_IgFLNa3_4m9p.dat", sep = "\t")
rog_l656f <- read.table("Rg_IgFLNa3_L656F.dat", sep = "\t")
rog_r484q <- read.table("Rg_IgFLNa3_R484Q.dat", sep = "\t")

rog_all <- cbind(rog_wt, rog_l656f[2], rog_r484q[2])
rog_all[1] <- seq(0, 99.99, 0.02)
colnames(rog_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
rog_seperate <- melt(rog_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rog_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(19.5, 22.5, 0.25))
p <- p + labs(x = "Time (ns)", y = "Rg (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("Rg_IgFLNa3_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Rg IgFLNa4 ---------------
rog_wt <- read.table("Rg_IgFLNa4_4m9p.dat", sep = "\t")
rog_l656f <- read.table("Rg_IgFLNa4_L656F.dat", sep = "\t")
rog_r484q <- read.table("Rg_IgFLNa4_R484Q.dat", sep = "\t")

rog_all <- cbind(rog_wt, rog_l656f[2], rog_r484q[2])
rog_all[1] <- seq(0, 99.99, 0.02)
colnames(rog_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
rog_seperate <- melt(rog_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rog_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(19.5, 22.5, 0.25))
p <- p + labs(x = "Time (ns)", y = "Rg (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("Rg_IgFLNa4_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Rg IgFLNa5 ---------------
rog_wt <- read.table("Rg_IgFLNa5_4m9p.dat", sep = "\t")
rog_l656f <- read.table("Rg_IgFLNa5_L656F.dat", sep = "\t")
rog_r484q <- read.table("Rg_IgFLNa5_R484Q.dat", sep = "\t")

rog_all <- cbind(rog_wt, rog_l656f[2], rog_r484q[2])
rog_all[1] <- seq(0, 99.99, 0.02)
colnames(rog_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
rog_seperate <- melt(rog_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rog_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(19.5, 22.5, 0.25))
p <- p + labs(x = "Time (ns)", y = "Rg (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("Rg_IgFLNa5_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Number of H-bonds btw IgFLNa3-4 ------------------
hbonds_wt <- read.table("hbondsBtwIgFLNa3_4_4m9p.dat", sep = " ")
hbonds_l656f <- read.table("hbondsBtwIgFLNa3_4_L656F.dat", sep = " ")
hbonds_r484q <- read.table("hbondsBtwIgFLNa3_4_R484Q.dat", sep = " ")

hbonds_all <- cbind(hbonds_wt, hbonds_l656f[2], hbonds_r484q[2])
hbonds_all[1] <- seq(0, 99.99, 0.02)
colnames(hbonds_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
hbonds_seperate <- melt(hbonds_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(hbonds_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 122, 10))
p <- p + labs(x = "Time (ns)", y = "Number of H-bonds")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("Hbonds_btwIgFLNa3-4_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Single residue RMSF for Trp582 -----------------
srRMSF_wt <- read.table("srRMSF_Trp582_4m9p.xvg", sep = "\t", skip = 8)[1:2]
srRMSF_l656f <- read.table("srRMSF_Trp582_L656F.xvg", sep = "\t", skip = 8)[1:2]
srRMSF_r484q <- read.table("srRMSF_Trp582_R484Q.xvg", sep = "\t", skip = 8)[1:2]

srRMSF_all <- cbind(srRMSF_wt, srRMSF_l656f[2], srRMSF_r484q[2])
srRMSF_all[1] <- seq(0, 99.99, 0.02)
colnames(srRMSF_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
srRMSF_seperate <- melt(srRMSF_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(srRMSF_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 22.5, 0.5), limits = c(0.5, 2.5))
p <- p + labs(x = "Time (ns)", y = "Single Residue RMSF (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.41, .99),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("srRMSF_Trp582_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Distance btw IgFLNa3-4 --------------------
dist_wt <- read.table("dist_IgFLNa3_4_4m9p.dat", sep = ",")
dist_l656f <- read.table("dist_IgFLNa3_4_L656F.dat", sep = ",")
dist_r484q <- read.table("dist_IgFLNa3_4_R484Q.dat", sep = ",")

dist_all <- cbind(dist_wt, dist_l656f[2], dist_r484q[2])
dist_all[1] <- seq(0, 99.99, 0.02)
colnames(dist_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
dist_seperate <- melt(dist_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(dist_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 22.5, 0.5))
p <- p + labs(x = "Time (ns)", y = "Distance (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("dist_IgFLNa3-4_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Distance btw IgFLNa3-5 --------------------
dist_wt <- read.table("dist_IgFLNa3_5_4m9p.dat", sep = ",")
dist_l656f <- read.table("dist_IgFLNa3_5_L656F.dat", sep = ",")
dist_r484q <- read.table("dist_IgFLNa3_5_R484Q.dat", sep = ",")

dist_all <- cbind(dist_wt, dist_l656f[2], dist_r484q[2])
dist_all[1] <- seq(0, 99.99, 0.02)
colnames(dist_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
dist_seperate <- melt(dist_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(dist_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 22.5, 0.5))
p <- p + labs(x = "Time (ns)", y = "Distance (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("dist_IgFLNa3-5_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Distance btw IgFLNa4-5 --------------------
dist_wt <- read.table("dist_IgFLNa4_5_4m9p.dat", sep = ",")
dist_l656f <- read.table("dist_IgFLNa4_5_L656F.dat", sep = ",")
dist_r484q <- read.table("dist_IgFLNa4_5_R484Q.dat", sep = ",")

dist_all <- cbind(dist_wt, dist_l656f[2], dist_r484q[2])
dist_all[1] <- seq(0, 99.99, 0.02)
colnames(dist_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
dist_seperate <- melt(dist_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(dist_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 22.5, 0.5))
p <- p + labs(x = "Time (ns)", y = "Distance (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("dist_IgFLNa4-5_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()

# Distance btw IgFLNa3-4 (only interaction surface residues) --------------------
dist_wt <- read.table("dist_IgFLNa3_4_onlyInteractionSurfaceResidues_4m9p.dat", sep = ",")
dist_l656f <- read.table("dist_IgFLNa3_4_onlyInteractionSurfaceResidues_L656F.dat", sep = ",")
dist_r484q <- read.table("dist_IgFLNa3_4_onlyInteractionSurfaceResidues_R484Q.dat", sep = ",")

dist_all <- cbind(dist_wt, dist_l656f[2], dist_r484q[2])
dist_all[1] <- seq(0, 99.99, 0.02)
colnames(dist_all) <- c("time","Wild-type", "p.Leu656Phe", "p.Arg484Gln")

library(reshape2)
dist_seperate <- melt(dist_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(dist_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 0.8) 
p <- p + geom_line(size=0.7)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 22.5, 0.5))
p <- p + labs(x = "Time (ns)", y = "Distance (Å)")
p <- p + ggsci::scale_color_jama()
p <- p + theme_minimal()
p <- p + theme(axis.title.y = element_text(angle = 90, size = 24),
               axis.title.x = element_text(size = 24),
               axis.ticks = element_line(),
               axis.text = element_text(size = 19),
               legend.title = element_blank(),
               legend.text = element_text(size = 15),
               legend.position = c(.01, .21),
               legend.justification = c("left", "top"),
               legend.box.just = "left",
               legend.margin = margin(6, 6, 6, 6),
               panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p

pdf("dist_IgFLNa3-4_onlyInteractionSurfaceResidues_4m9pCases.pdf", width = 8, height = 6)
p
dev.off()