#!/usr/bin/Rscript
rm(list=ls())
library(ggplot2)
alldata =c()
for(i in 1:10){
filename=paste("../Fig3/data/H_cycle16_LG", i, ".txt", sep="")
Hdata <- read.table(filename, header=T)
alldata <- rbind(alldata, Hdata)
}

filename="../data/2cm_sim_collated_results.txt"
simdata <- read.table(filename, header=T)


sigs=simdata[!is.na(simdata$Sigat001),]

sigs$Sigat001 = factor(sigs$Sigat001, ordered=T, levels=c("BSCB", "BSSS", "Both"))
levels(sigs$Sigat001) = c("BSCB1", "BSSS", "Both")
blocks=geom_rect(data=sigs, aes(fill=Sigat001, xmin=phystart/1000000, xmax=phystop/1000000, ymin=0, ymax=0.5)) 
BSSS_H <- geom_point(data=alldata, aes(x=physpos/1000000, y=HetSS), colour="blue", size=1)
BSCB1_H <- geom_point(data=alldata, aes(x=physpos/1000000, y=HetNSS), colour="red2", size=1)
png('../FigS2A.png')
ggplot(data=alldata) + blocks + BSSS_H + BSCB1_H + facet_grid(LG~.) + ylab("Heterozygosity") + xlab("Position (Mb)") + theme_bw() +
scale_fill_manual(values=c("pink1", "skyblue2", "green"), name="P < 0.001") + theme(panel.border=element_rect(size=0.1))  + theme(axis.text.y=element_text(size=8))
dev.off()

blocks=geom_rect(data=sigs, aes(fill=Sigat001, xmin=genstart, xmax=genstop, ymin=0, ymax=0.5)) 
BSSS_H <- geom_point(data=alldata, aes(x=genpos, y=HetSS), colour="blue", size=1)
BSCB1_H <- geom_point(data=alldata, aes(x=genpos, y=HetNSS), colour="red2", size=1)
png('../FigS2B.png')
ggplot(data=alldata) + blocks + BSSS_H + BSCB1_H + facet_grid(LG~.) + ylab("Heterozygosity") + xlab("Position (cM)") + theme_bw() +
  scale_fill_manual(values=c("pink1", "skyblue2", "green"), name="P < 0.001") + theme(panel.border=element_rect(size=0.1))  + theme(axis.text.y=element_text(size=8))
dev.off()

