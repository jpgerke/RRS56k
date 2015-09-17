#!/usr/bin/Rscript
rm(list=ls())
library(ggplot2)
library(dplyr)
alldata =c()
for(i in 1:10){
filename=paste("./data/H_cycle16_LG", i, ".txt", sep="")
Hdata <- read.table(filename, header=T)
alldata <- rbind(alldata, Hdata)
}

BSSSHdata = alldata[,-2]
BSCBHdata = alldata[,-1]
colnames(BSSSHdata) = c("H", "LG", "Cycle", "Start", "Stop", "physpos", "Genpos", "Rate")
colnames(BSCBHdata) = c("H", "LG", "Cycle", "Start", "Stop", "physpos", "Genpos", "Rate")
BSSSHdata$Population = "BSSS"
BSCBHdata$Population = "BSCB1"
both = rbind(BSSSHdata, BSCBHdata)

png("../fig3.png", height=900, width=700)
ggplot(data=both) + aes(x=physpos/1000000, y=H, colour=Population) + facet_grid(LG~Population) + geom_point(size=1) + scale_colour_manual(values=c("blue", "red2")) + theme_bw() + xlab("Physical Position (MB)")
dev.off()


jpeg("../Fig_S2.jpg", quality=1000,  width=700, height=900)
ggplot(data=both) + aes(x=physpos/1000000, y=H, colour=Population) + facet_grid(LG~Population) + geom_point(size=1) + scale_colour_manual(values=c("blue", "red2")) + theme_bw() + xlab("Physical Position (MB)")
dev.off()
