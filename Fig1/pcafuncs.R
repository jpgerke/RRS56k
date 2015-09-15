rm(list=ls())
load("../data/new_full_set.RData")
source("../rcode/writing_functions.R")
source("../rcode/building_functions.R")
source("../rcode/formatting_functions.R")


convert_numeric = function(allele){
  calls = unique(allele)[!is.na(unique(allele))]
  is.homozygote = sapply(calls, function(x){
    myalleles = unlist(strsplit(x,''))
    boolean = myalleles[1] == myalleles[2]
    return(boolean)
  })
  allele[allele == calls[is.homozygote][1]] = 0 
  allele[allele == calls[is.homozygote][2]] = 2
  allele[allele == calls[!is.homozygote][1]] = 1
  return(as.numeric(allele))
  }
  

mynumeric = t(apply(new_full_set$calls,1,convert_numeric))
rownames(mynumeric) = rownames(new_full_set$calls)
colnames(mynumeric) = colnames(new_full_set$calls)

prune_monos = function(datamat){
  polymorphics = apply(datamat, 1, function(x){sum(!is.na(unique(x)))})
  keeps = polymorphics>1
  newmat = datamat[keeps,]
  return(newmat)
}

centermean = function(datamat){
  mymeans = getmeans(datamat)
  centered = (datamat - mymeans)
  centered[is.na(centered)] = 0
  return(list(centered, mymeans))
}

getmeans = function(datamat){
  means = apply(datamat, 1, mean, na.rm=T)
  return(means)
}

runPCA = function(datamat){
  centerfunc = centermean(datamat)
  means = centerfunc[[2]]
  centered = centerfunc[[1]]
  centered = t(centered)
  centered = centered/(sqrt(ncol(centered)))
  PCs = prcomp(centered, center=FALSE, retx=T)
  loadings = t(t(PCs$rotation) *PCs$sdev)
  projection = centered %*% loadings
  varexp = (PCs$sdev^2)/sum(PCs$sdev^2)
  results = list(M=centered, projections=projection, loadings=loadings, eigs = PCs$sdev^2, means = means, nmarks=ncol(centered), varexp = varexp)
  return(results)
}

make_projection = function(PCAobj,newdat){
  newdat = t(newdat)
  tokeep = colnames(newdat) %in% colnames(PCAobj$M)
  pruned = newdat[,tokeep]
  centered = t(t(pruned) - PCAobj$means)
  centered[is.na(centered)] = 0
  centered = centered/(sqrt(PCAobj$nmarks))
  projection = centered %*% PCAobj$loadings
  return(list(M=centered, projection=projection))
}

test = mynumeric[1:10,1:5]
test = prune_monos(test)
mine = runPCA(test)
mine2 = make_projection(mine, test)
library(ggplot2)

PCAset <- remove.line(line="Ill_Hy_SS",dataset=new_full_set)
outbreds = PCAset$calls[,PCAset$status=='Outbred']
outbred_num =t(apply(outbreds,1,convert_numeric))
outbred_num = prune_monos(outbred_num)

out = runPCA(outbred_num)

founds = PCAset$calls[,PCAset$status=='Founder']
found_num = t(apply(founds, 1, convert_numeric))
proj = make_projection(out, found_num)

plot(proj$projection[,1], proj$projection[,2])

rownames(out$projections) = colnames(outbreds)
rownames(proj$projection) = colnames(founds)

outdf = data.frame(Inbred=rownames(out$projections), pc1 = out$projections[,1], pc2 = out$projections[,2], Group=new_full_set$cycle[match(rownames(out$projections), colnames(new_full_set$calls))], hetgrp=new_full_set$group[match(rownames(out$projections), colnames(new_full_set$calls))] )
projdf = data.frame(Inbred=rownames(proj$projection), pc1 = proj$projection[,1], pc2 = proj$projection[,2], Group=new_full_set$cycle[match(rownames(proj$projection), colnames(new_full_set$calls))], hetgrp=new_full_set$group[match(rownames(proj$projection), colnames(new_full_set$calls))])

projdf$Grouping = ''
projdf[projdf$Group =='Founder' & projdf$hetgrp=="BSSS",]$Grouping = 'BSSS Founder'
projdf[projdf$Group =='Founder' & projdf$hetgrp=="BSCB",]$Grouping = 'BSCB1 Founder'
outdf$Grouping = as.character(outdf$Group)
projdf$Grouping = as.character(projdf$Grouping)

fulldf = rbind(outdf,projdf)
fulldf$Cycle = factor(fulldf$Grouping, levels=c("BSSS Founder", "BSCB1 Founder", "0", "4", "8", "12", "16"), ordered=T)

library(ggplot2)
library(RColorBrewer)
library(cowplot)
mycols = brewer.pal(8, "Set1")


testplot = ggplot(data=fulldf) + aes(x=pc1, y=-pc2, colour=Cycle) + geom_point(size=2, alpha=0.8)  + scale_colour_manual(values=mycols)
testplot= testplot + theme_bw() + ylab("PC2 (7%)") + xlab("PC1 (30%)") + theme(axis.title.x=element_text(vjust=-0.5,hjust=0.6))

jpeg('fig1.jpg', quality = 400, width=6, height=4, units="in", res=500)
testplot + geom_text(aes(x=x1, y=x2, label=texthere), data.frame(x1=c(-.13,.12), x2=c(0.03,0.03), texthere=c("BSCB1", "BSSS")),colour="black")
dev.off()

ggplot(data=outdf) + aes(x=pc1, y=-pc2, colour=as.factor(hetgrp)) + geom_point() + geom_point(data=projdf, aes(label=Inbred))
