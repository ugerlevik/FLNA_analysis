setwd("D:/SezermanLab/FLNA/cMDs/data/1_4m9p_cases/")
library(ggplot2)
library(ggpubr)

palet <- ggsci::pal_startrek()
colorPalet <- palet(3)

# hbonds (focused) -------------
hbonds_wt <- read.table("hbonds10ofresid 484_wt.dat", sep = " ")
hbonds_r484q <- read.table("hbonds10ofresid 484_R484Q.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_r484q[2])
hbonds_all[1] <- seq(0, 99.99, 0.02)
colnames(hbonds_all) <-  c("time","WT","Arg484Gln")
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p1 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p1 <- p1 + geom_boxplot(size = 1.1, alpha = 0.75)
# p1 <- p1 + geom_line(size = 1, alpha = 0.75)
# p1 <- p1 + scale_x_continuous(breaks = seq(0, 100, 10))
p1 <- p1 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,25))
p1 <- p1 + labs(x = "Time (ns)", y = "Number of H-bonds")
p1 <- p1 + scale_color_manual(values = colorPalet[c(1, 2)])
p1 <- p1 + theme_minimal()
p1 <- p1 + theme(#axis.title.y = element_text(angle = 90, size = 14),
  # axis.title.x = element_text(size = 14),
  axis.title = element_blank(),
  axis.ticks = element_line(),
  axis.text = element_text(size = 16),
  axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 14),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = "none",
  legend.key.width = unit(1.25,"cm"),
  panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p1

hbonds_wt <- read.table("hbonds10ofresid 656_wt.dat", sep = " ")
hbonds_l656f <- read.table("hbonds10ofresid 656_L656F.dat", sep = " ")

hbonds_all <- data.frame(hbonds_wt, hbonds_l656f[2])
hbonds_all[1] <- seq(0, 99.99, 0.02)
colnames(hbonds_all) <-  c("time","WT","Leu656Phe")
hbonds_seperate <- reshape2::melt(hbonds_all[1:3], id.var = "time") 

p2 <- ggplot(hbonds_seperate, aes(x=variable, y=value, colour=variable))
p2 <- p2 + geom_boxplot(size = 1.1, alpha = 0.75)
# p2 <- p2 + geom_line(size = 1, alpha = 0.75)
# p2 <- p2 + scale_x_continuous(breaks = seq(0, 100, 10))
p2 <- p2 + scale_y_continuous(breaks = seq(0, 100, 2), limits = c(0,25))
p2 <- p2 + labs(x = "Time (ns)", y = "Number of H-bonds")
p2 <- p2 + scale_color_manual(values = colorPalet[c(1, 3)])
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(#axis.title.y = element_text(angle = 90, size = 14),
  # axis.title.x = element_text(size = 14),
  axis.title = element_blank(),
  axis.ticks.x = element_line(),
  axis.ticks.y = element_line(),
  axis.text.y = element_blank(),
  axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 14),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.position = "none",
  legend.key.width = unit(1.25,"cm"),
  panel.border = element_rect(colour = "gray60", fill = NA, size = .5))
p2

fig <- ggarrange(p1,p2, widths = c(2,2),
                 nrow = 1, ncol = 2, align = "h")

pdf("hbondsWithin10_4m9pCases.pdf", width = 4, height = 3)
annotate_figure(fig, left = text_grob("Number of H-bonds", size = 14, rot = 90))
dev.off()

