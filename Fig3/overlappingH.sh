#!/bin/bash

for cycle in 0 4 8 12 16
do

for chrom in {1..10}
do

R --slave --no-save <<HEREFILE

source("../rcode/building_functions.R")
source("../rcode/analysis_functions.R")
source("../rcode/formatting_functions.R")
source("../rcode/writing_functions.R")
load("../data/new_full_set.RData")



	currentseq=1:sum(new_full_set\$LG==${chrom})
        windows <- define.overlapping.windows(currentseq, window=15, step=5)

	chrommarkers <- new_full_set\$calls[new_full_set\$LG==${chrom},]
	chromgen <- new_full_set\$F2map[new_full_set\$LG==${chrom}]
	chromphys <- new_full_set\$position[new_full_set\$LG==${chrom}]

windowtable <- t(apply(windows, 1, function(x){
  midmarker <- ceiling(median(x[1]:x[2]))
  physpos <- chromphys[midmarker]
  phystart <- chromphys[x[1]]
  phystop <- chromphys[x[2]]
  genstart <- chromgen[x[1]]
  genstop <- chromgen[x[2]]
  genpos <- chromgen[midmarker]
  rate <-(genstop-genstart)/(phystop/1000000 - phystart/1000000)
  return(c(x[1], x[2], physpos, genpos, rate))  
  }))

  windowtable <- data.frame(windowtable)

  colnames(windowtable) <- c("start", "stop", "physpos", "genpos", "rate")

  outbreds <- split.status(new_full_set, "Outbred")

outbreds <- split.all(outbreds)

realresults <- t(apply(windows, 1, function(x){
  BSSScalls <- outbreds\$BSSS_${cycle}\$calls[outbreds\$BSSS_${cycle}\$LG==${chrom},]  
  BSCBcalls <- outbreds\$BSCB_${cycle}\$calls[outbreds\$BSCB_${cycle}\$LG==${chrom},]
  hetSS <- averageH(BSSScalls[x[1]:x[2],])
  hetNSS <- averageH(BSCBcalls[x[1]:x[2],])
  return(c(hetSS, hetNSS, ${chrom}, ${cycle}))
}))

colnames(realresults) <- c("HetSS", "HetNSS", "LG", "Cycle")
realresults <- cbind(realresults, windowtable)


write.table(realresults, file="./data/H_cycle${cycle}_LG${chrom}.txt", sep="\t", quote=F, row.names=F, col.names=T)

HEREFILE
done
done
