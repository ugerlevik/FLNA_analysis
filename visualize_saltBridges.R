#########################################################################
## Project: FLNA                                 
## Script purpose: Visualize salt bridges
## Date: Mar 30, 2020                                                   
## Author: Umut Gerlevik                                                 
#########################################################################
options(stringsAsFactors = FALSE)
library(ggplot2)
library(reshape2)

# find the residues:
# set a [atomselect top "charged and name CA and same residue as (within 10 of (protein and resid 127))"]
# $a get {resname resid chain}

# dcdlist <- c("4m9p", "L656F")
# saltbrList <- c()
# for (i in dcdlist) {
#   saltbrList <- c(saltbrList, dir(paste("data/3_3hop_cases/saltbridges_", i, sep = ""), pattern = ".dat"))
# }
# saltbrList <- unique(saltbrList)

# w10_4m9p_484 <- c("ARG484", "LYS504", "GLU642", "ARG664",
#                   "ASP665", "ARG563", "ASP535", "LYS508",
#                   "GLU639", "GLU587", "ASP662")
# w10_4m9p_656 <- c("LYS498", "ARG655", "LYS493", "GLU499",
#                   "GLU587", "GLU715", "ARG488", "ARG580",
#                   "ARG496", "ASP502", "GLU652", "ASP653")
# w10_3cnk_2593 <- c("LYS2631_chainA", "LYS2631_chainB", "GLU2602_chainB", "ARG2643_chainB",
#                    "GLU2635_chainA", "ARG2643_chainA", "GLU2603_chainB")
# GLU2593_saltbr <- dir("data/3_3hop_cases/saltbridges_L656F/", pattern = ".dat")
# GLU2593_saltbr <- GLU2593_saltbr[grep("GLU2593", GLU2593_saltbr)] # mutated to GLU
# w10_3cnk_2641 <- c("GLU2635_chainA", "GLU2635_chainB", "ARG2643_chainA", "ARG2643_chainB",
#                    "LYS2631_chainB", "GLU2625_chainA", "LYS2631_chainA")
w10_3cnk_2622_2623 <- c("LYS2569_chainA", "LYS2575_chainA", "ARG2598_chainA", "GLU2602_chainA",
                        "GLU2603_chainA", "ARG2612_chainA", "LYS2621_chainA", "ASP2622_chainA",
                        "LYS2623_chainA", "GLU2625_chainA", "ARG2598_chainB", "GLU2602_chainB",
                        "GLU2603_chainB", "ARG2612_chainB", "LYS2621_chainB", "ASP2622_chainB",
                        "LYS2623_chainB", "GLU2625_chainB", "ARG2643_chainB")
# w10_3hop_39 <- c("GLU37", "ASP38", "LYS42", "LYS135", "ASP249")
# GLU39_saltbr <- dir("data/3_3hop_cases/saltbridges_L656F/", pattern = ".dat")
# GLU39_saltbr <- GLU39_saltbr[grep("GLU39", GLU39_saltbr)] # mutated to GLU
# w10_3hop_119 <- c("LYS87", "GLU112", "ASP115", "ARG116", "GLU117",
#                   "LYS120", "ASP156", "GLU157")
# w10_3hop_80 <- c("ARG51", "GLU55", "LYS58", "LYS62", "ARG63",
#                  "ASP70", "ASP73", "ARG76", "GLU82", "LYS87",
#                  "LYS88")
# w10_3hop_102 <- c("ASP73", "ARG91", "LYS92", "ARG96",
#                   "ARG100","GLU105", "LYS127")
# w10_3hop_149 <- c("GLU55", "LYS169", "ARG171", "GLU254", "LYS266")
# w10_3hop_129 <- c("ASP70", "ASP73", "ARG76", "ARG100", "GLU105",
#                   "GLU112", "ASP125", "LYS127", "ASP131", "LYS135")
# w10_3hop_127 <- c("ASP73","ARG100","GLU105","ASP125","LYS127","ASP131")
# w10_2brq_2257 <- c("GLU2249_chainA","ARG2250_chainA","GLU2252_chainA",
#                    "GLU2258_chainA","LYS2280_chainA","GLU2282_chainA",
#                    "GLU2286_chainA","ARG2288_chainA","GLU2301_chainA",
#                    "GLU2249_chainB","ARG2250_chainB","GLU2252_chainB",
#                    "GLU2258_chainB","GLU2276_chainB","GLU2282_chainB",
#                    "GLU2286_chainB","ARG2288_chainB","GLU2301_chainB",
#                    "ASP2304_chainB","GLU2306_chainB")

saltbrFileName <- c()
for (i in w10_3cnk_2622_2623) {
  saltbrFileName <- c(saltbrFileName, dir("data/2_3cnk_cases/saltbridges_3cnk/", pattern = ".dat")[grep(i, dir("data/2_3cnk_cases/saltbridges_3cnk/", pattern = ".dat"))])
  saltbrFileName <- c(saltbrFileName, dir("data/2_3cnk_cases/saltbridges_D2622_K2623del/", pattern = ".dat")[grep(i, dir("data/2_3cnk_cases/saltbridges_D2622_K2623del/", pattern = ".dat"))])
}
# saltbrFileName <- c(saltbrFileName, GLU2593_saltbr)
# saltbrFileName <- c(saltbrFileName, GLU39_saltbr)

for (j in saltbrFileName) {
  if (j %in% dir("data/2_3cnk_cases/saltbridges_3cnk/") &&
      j %in% dir("data/2_3cnk_cases/saltbridges_D2622_K2623del/")) {
    
    wt <- read.table(paste("data/2_3cnk_cases/saltbridges_3cnk/", j, sep = ""), sep = " ")
    mut <- read.table(paste("data/2_3cnk_cases/saltbridges_D2622_K2623del/", j, sep = ""), sep = " ")
    wt_mut <- cbind(wt, mut[2])
    wt_mut[1] <- seq(0, 99.99, 0.02)
    colnames(wt_mut) <- c("Time", "Wild-type", "D2622_K2623del")
    wt_mut <- melt(wt_mut[1:3], id.var = "Time")
    p <- ggplot(wt_mut, aes(x = Time, y = value, colour = variable))
    p <- p + geom_point(size = 1)
    p <- p + geom_line(size = 0.8)
    p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
    p <- p + scale_y_continuous(breaks = seq(0, 30, 1))
    p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
    p <- p + ggsci::scale_color_jama()
    p <- p + theme_minimal()
    p <- p + theme(axis.title.y = element_text(angle = 90, size = 20),
                   axis.title.x = element_text(size = 24),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 19),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = c(.01, 1.01),
                   legend.justification = c("left", "top"),
                   legend.box.just = "left",
                   legend.margin = margin(6, 6, 6, 6),
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
    pdf(file = paste("output/saltbr/D2622_K2623del/WT_D2622_K2623del_", sub(".dat", "", j), ".pdf", sep = ""), width = 8, height = 6)
    print(p)
    dev.off()
    
  } else if ((j %in% dir("data/2_3cnk_cases/saltbridges_3cnk/")) &&
             !(j %in% dir("data/2_3cnk_cases/saltbridges_D2622_K2623del/"))) {
    
    wt <- read.table(paste("data/2_3cnk_cases/saltbridges_3cnk/", j, sep = ""), sep = " ")
    wt[1] <- seq(0, 99.99, 0.02)
    colnames(wt) <- c("Time", "saltbr")
    p <- ggplot(wt, aes(x = Time, y = saltbr))
    p <- p + geom_point(size = 1)
    p <- p + geom_line(size = 0.8)
    p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
    p <- p + scale_y_continuous(breaks = seq(0, 30, 1))
    p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
    p <- p + ggsci::scale_color_jama()
    p <- p + theme_minimal()
    p <- p + theme(axis.title.y = element_text(angle = 90, size = 20),
                   axis.title.x = element_text(size = 24),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 19),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = c(.01, 1.01),
                   legend.justification = c("left", "top"),
                   legend.box.just = "left",
                   legend.margin = margin(6, 6, 6, 6),
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
    pdf(file = paste("output/saltbr/D2622_K2623del/WT_", sub(".dat", "", j), ".pdf", sep = ""), width = 8, height = 6)
    print(p)
    dev.off()
    
  } else if (!(j %in% dir("data/2_3cnk_cases/saltbridges_3cnk/")) &&
             j %in% dir("data/2_3cnk_cases/saltbridges_D2622_K2623del/")) {
    
    mut <- read.table(paste("data/2_3cnk_cases/saltbridges_D2622_K2623del/", j, sep = ""), sep = " ")
    mut[1] <- seq(0, 99.99, 0.02)
    colnames(mut) <- c("Time", "saltbr")
    p <- ggplot(mut, aes(x = Time, y = saltbr))
    p <- p + geom_point(size = 1)
    p <- p + geom_line(size = 0.8)
    p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
    p <- p + scale_y_continuous(breaks = seq(0, 30, 1))
    p <- p + labs(x = "Time (ns)", y = paste(sub("saltbr-", "", sub(".dat", "", j)), "(Å)", sep = " "))
    p <- p + ggsci::scale_color_jama()
    p <- p + theme_minimal()
    p <- p + theme(axis.title.y = element_text(angle = 90, size = 20),
                   axis.title.x = element_text(size = 24),
                   axis.ticks = element_line(),
                   axis.text = element_text(size = 19),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 15),
                   legend.position = c(.01, 1.01),
                   legend.justification = c("left", "top"),
                   legend.box.just = "left",
                   legend.margin = margin(6, 6, 6, 6),
                   panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
    pdf(file = paste("output/saltbr/D2622_K2623del/D2622_K2623del_", sub(".dat", "", j), ".pdf", sep = ""), width = 8, height = 6)
    print(p)
    dev.off()
    
  }
}
rm(list = ls())
