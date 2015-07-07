#!/usr/bin/Rscript

rm(list=ls())
library(ggplot2)
library(dplyr)

loadcycle = function(cycle, mychrom){
  infile = paste("../Fig3/data/H_cycle", cycle, "_LG", mychrom, ".txt", sep='')
  myfile = read.table(infile, header=T)
  return(myfile)
}
  
loadchrom = function(chrom){
  cycles = c(0,4,8,12,16)
  tables = lapply(cycles, loadcycle, mychrom=chrom) 
  results = do.call(rbind, tables)
  return(results)
}

chromdata_H = lapply(1:10, loadchrom)

simresults = read.table("../data/2cm_sim_collated_results.txt", header=T)


#The code to make plots is ugly because ggplot has some unpredictable (at least to me) behavior when wrapped in functions due to scope issues.
####code to make plot for chromosome 9 BSSS.  ####

  Chrom=9 
  Hetgrp="BSSS"
  Hdata = chromdata_H[[Chrom]]
  mytitle = paste("Chromosome ", Chrom, ", ", Hetgrp, sep="" )
  
  sigregions001 = filter(simresults, Sigat001 == "BSSS" | Sigat001 == "Both") %>% filter(LG==Chrom)
  sigregions0025 = filter(simresults, Sigat0025 == "BSSS" | Sigat0025 == "Both") %>% filter(LG==Chrom)
	sigat001rect_phys <-  geom_rect(data=sigregions001, aes(xmin=phystart/1000000, xmax=phystop/1000000, ymin=0, ymax=0.5), fill="lightskyblue") 
#	sigat0025rect_phys <- geom_rect(data=sigregions0025, aes(xmin=phystart/1000000, xmax=phystop/1000000, ymin=0, ymax=0.5), fill="lightskyblue") 
	Hplot_phys <- geom_point(size = 1, data=Hdata, aes(x=physpos/1000000, y=HetSS, colour=sqrt(rate^(1/4)))) 
	myphysplot=ggplot(data=Hdata) + facet_grid(Cycle~.) + sigat001rect_phys + Hplot_phys  + xlab("Position(Mb)")+ ylab("Heterozygosity") + scale_colour_gradient2(low="red", mid="orange", high="yellow", midpoint=median(Hdata$rate^(1/4))) + theme_bw() + ggtitle("BSSS, Chromosome 9")

  axesopts = theme(axis.text.y=element_text(size=6), axis.text.x=element_text(size=6), axis.title.y=element_text(size=8, angle=90), axis.title.x=element_text(size=8))

  myphysplot = myphysplot + axesopts + theme(legend.position="None")  #+ ggtitle(mytitle)

  simdata1rect_gen <-  geom_rect(data=sigregions001, aes(xmin=genstart, xmax=genstop, ymin=0, ymax=0.5), fill="lightskyblue") 
 # simdata2rect_gen <- geom_rect(data=sigregions0025, aes(xmin=genstart, xmax=genstop, ymin=0, ymax=0.5), fill="lightskyblue") 
  Hplot_gen <- geom_point(size = 1, data=Hdata, aes(x=genpos, y=HetSS, colour=sqrt(rate^(1/4)))) 
  mygenplot=ggplot(data=Hdata) + facet_grid(Cycle~.) + simdata1rect_gen  + Hplot_gen + xlab("Position(cM)")+ ylab("Heterozygosity") + scale_colour_gradient2(low="red", mid="orange", high="yellow", midpoint=median(Hdata$rate^(1/4))) + theme_bw() #+ ggtitle("BSSS Chromosome 9")
  
  mygenplot = mygenplot + axesopts + theme(legend.position="None")
  
  Chr9BSSS = list(phys=myphysplot, gen=mygenplot)

 

# outfile_phys="../Fig4A.pdf"
# pdf(outfile_phys, width=5,height=4)
# print(Chr9BSSS[["phys"]])
# dev.off()
# 
# outfile="../Fig4B.pdf"
# pdf(outfile, width=5,height=4)
# print(Chr9BSSS[["gen"]])
# dev.off()


#### now chr4 BSCB1
Chrom=4 
Hetgrp="BSCB1"
Hdata = chromdata_H[[Chrom]]
mytitle = paste("Chromosome ", Chrom, ", ", Hetgrp, sep="" )

sigregions001 = filter(simresults, Sigat001 == "BSCB" | Sigat001 == "Both") %>% filter(LG==Chrom)
sigregions0025 = filter(simresults, Sigat0025 == "BSCB" | Sigat0025 == "Both") %>% filter(LG==Chrom)
sigat001rect_phys <-  geom_rect(data=sigregions001, aes(xmin=phystart/1000000, xmax=phystop/1000000, ymin=0, ymax=0.5), fill="lightskyblue") 
#sigat0025rect_phys <- geom_rect(data=sigregions0025, aes(xmin=phystart/1000000, xmax=phystop/1000000, ymin=0, ymax=0.5), fill="lightskyblue") 
Hplot_phys <- geom_point(size = 1, data=Hdata, aes(x=physpos/1000000, y=HetNSS, colour=sqrt(rate^(1/4)))) 
myphysplot=ggplot(data=Hdata) + facet_grid(Cycle~.) + sigat001rect_phys + Hplot_phys + xlab("Position(Mb)")+ ylab("Heterozygosity") + scale_colour_gradient2(low="red", mid="orange", high="yellow", midpoint=median(Hdata$rate^(1/4))) + theme_bw() + ggtitle("BSCB1, Chromosome 4")

axesopts = theme(axis.text.y=element_text(size=6), axis.text.x=element_text(size=6), axis.title.y=element_text(size=8, angle=90), axis.title.x=element_text(size=8))

myphysplot = myphysplot + axesopts + theme(legend.position="None")  #+ ggtitle(mytitle)

simdata1rect_gen <-  geom_rect(data=sigregions001, aes(xmin=genstart, xmax=genstop, ymin=0, ymax=0.5), fill="lightskyblue") 
#simdata2rect_gen <- geom_rect(data=sigregions0025, aes(xmin=genstart, xmax=genstop, ymin=0, ymax=0.5), fill="lightskyblue") 
Hplot_gen <- geom_point(size = 1, data=Hdata, aes(x=genpos, y=HetNSS, colour=sqrt(rate^(1/4)))) 
mygenplot=ggplot(data=Hdata) + facet_grid(Cycle~.) + simdata1rect_gen + Hplot_gen  + xlab("Position(cM)")+ ylab("Heterozygosity") + scale_colour_gradient2(low="red", mid="orange", high="yellow", midpoint=median(Hdata$rate^(1/4))) + theme_bw()#+ ggtitle("BSCB1 Chromosome 4")n 

mygenplot = mygenplot + axesopts + theme(legend.position="None")

Chr4BSCB1 = list(phys=myphysplot, gen=mygenplot)

# outfile_phys="../Fig4C.pdf"
# pdf(outfile_phys, width=5,height=4)
# print(Chr4BSCB1[["phys"]])
# dev.off()
# 
# outfile="../Fig4D.pdf"
# pdf(outfile, width=5,height=4)
# print(Chr4BSCB1[["gen"]])
# dev.off()

library(cowplot)
allfour = plot_grid(Chr9BSSS[["phys"]], Chr4BSCB1[["phys"]], Chr9BSSS[["gen"]], Chr4BSCB1[["gen"]], labels=c("A", "B", "C", "D"), ncol=2)
save_plot("../fig4_combined.pdf", allfour,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)

