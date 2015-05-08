#!/usr/bin/Rscript
rm(list=ls())
library(ggplot2)
library(reshape)



simresults <- read.table("../data/BSSS_fullsim_HH.txt")
colnames(simresults) = c("Cycle4", "Cycle8", "Cycle12", "Cycle16")
melted_CB_all = melt(simresults)
colnames(melted_CB_all) = c("Cycle", "H")
melted_CB_all$Pruning = "all"

colnames(melted_CB_all) = c("Cycle",  "H", "Pruning")
medians = plyr::ddply(.data=melted_CB_all, .variables=c("Cycle", "Pruning"), .fun=function(x){median(x[,2])})
colnames(medians) = c("Cycle", "Pruning",  "Med")
realdata = read.table("BSSS_realresults.txt", header=T)
getbounds = function(x){
  quantile(x[,2], probs=c(0.01, 0.99))
}
bounds = plyr::ddply(.data=melted_CB_all, .variables=c("Cycle", "Pruning"), .fun=getbounds)
colnames(bounds) = c("Cycle", "Pruning", "one", "ninenine")
realdata = realdata[realdata$Pruning=="all",]
png("../FigS1A.png", width=400,height=300)
myplot = ggplot()
myplot = myplot + geom_point(data=melted_CB_all, aes(x=Cycle, y=H), alpha=0.1) + geom_line(data=bounds, aes(x=Cycle, y=one, group=Pruning), color="black") + geom_line(data=bounds, aes(x=Cycle, y=ninenine, group=Pruning), color="black")
myplot = myplot + geom_point(data=medians, aes(x=Cycle, y=Med), color="green") 
myplot =myplot + geom_line(data=realdata, aes(x=Cycle, y=H, group=Pruning), color="red") + geom_point(data=realdata, aes(x=Cycle, y=H), color="red") + theme_bw() + ggtitle("BSSS Genome-wide Simulation")
myplot + theme(axis.title.x=element_text(size=0))
dev.off()

# png("BSSS_sims_gray.png", width=1000,height=600)
# test = ggplot()
# test = test + geom_point(data=allresults, aes(x=Cycle, y=H), alpha=0.1) + geom_line(data=bounds, aes(x=Cycle, y=one, group=Pruning), color="black") + geom_line(data=bounds, aes(x=Cycle, y=ninenine, group=Pruning), color="black")+ facet_grid(.~Pruning)
# test = test + geom_point(data=medians, aes(x=Cycle, y=Med), color="green") 
# test + geom_line(data=realdata, aes(x=Cycle, y=H, group=Pruning), color="red") + geom_point(data=realdata, aes(x=Cycle, y=H), color="red") +facet_grid(.~Pruning) + opts(title="BSSS Genome-wide Sims")
# dev.off()
