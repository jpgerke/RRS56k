library(dplyr)

data = read.table("./data/2cM_sim_collated_results.txt", header=T)
mydata = select(data, Cycle, LG, phystart, phystop, 
                genstart, genstop, unph_H_BSSS, unph_H_BSCB, 
                H_BSSS, H_BSCB, BSSS_quantile_0.001, BSCB_quantile_0.001)
write.table(mydata, file="File_S2.txt", row.names=F, quote=F)
