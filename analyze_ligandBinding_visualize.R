setwd("D:/SezermanLab/FLNA/cMDs/data/1_4m9p_cases/")

# RMSD ------------
rmsd_wt <- read.table("rmsd_ligandBinding_4m9p.dat", sep = " ")
rmsd_l656f <- read.table("rmsd_ligandBinding_L656F.dat", sep = " ")
rmsd_r484q <- read.table("rmsd_ligandBinding_R484Q.dat", sep = " ")

rmsd_all <- cbind(rmsd_wt, rmsd_r484q[2], rmsd_l656f[2])
rmsd_all[1] <- seq(0, 99.99, 0.02)
colnames(rmsd_all) <- c("time","Wild-type", "p.Arg484Gln", "p.Leu656Phe")

rmsd_seperate <- reshape2::melt(rmsd_all[1:4], id.var = "time") 
library(ggplot2)
p <- ggplot(rmsd_seperate, aes(x=time, y=value, colour=variable)) 
p <- p + geom_point(size = 1.1, alpha = 0.75)
p <- p + geom_line(size = 1, alpha = 0.75)
p <- p + scale_x_continuous(breaks = seq(0, 100, 10))
p <- p + scale_y_continuous(breaks = seq(0, 7, 0.5), limits = c(0,2))
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

pdf("RMSD_ligandBinding_4m9pCases.pdf", width = 10, height = 6)
p
dev.off()
