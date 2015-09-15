rm(list=ls())

data = read.csv("Fig2_stats.csv")
data$Stat = c(rep("Fst", 6), rep("H", 18))
library(ggplot2)
data$Cycle = factor(data$Cycle,levels=c("Founders", "0", "4", "8", "12", "16"), ordered=TRUE)
data$Stat = factor(data$Stat,levels=c("H", "Fst"), ordered=TRUE)
data$Measure = as.character(data$Measure)
data$Measure[data$Measure == 'H BSCB'] = "H BSCB1"
jpeg("Fig2.jpg", height=5, width=7, units="in", res=500)
ggplot(data=data) + aes(x=Cycle, y=Value, group=Measure, colour=Measure) + geom_point(size=4) + facet_grid(.~Stat) + ylim(0,0.5) + geom_path(size=1) + theme_bw() + scale_colour_discrete("Statistic")
dev.off()
