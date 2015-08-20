rm(list = ls())

load("../data/new_full_set.RData")

sample_info <- data.frame(Name = new_full_set$lines,
                          Population = new_full_set$group,
                          Type = new_full_set$status,
                          Cycle = new_full_set$cycle)

sample_info$Type <- as.character(sample_info$Type)
sample_info$Type[sample_info$Type == "Outbred"] = "Population Sample"

markerinfo <- data.frame(Name = new_full_set$markernames,
                         ID = new_full_set$markercodes,
                         Chromosome = new_full_set$LG,
                         Position = new_full_set$position,
                         Genetic_Position = new_full_set$F2map)

write.csv(sample_info, file = "sample_info.csv", quote = F, row.names = F)
write.csv(markerinfo, file = "SNP_info.csv", quote = F, row.names = F)
write.table(new_full_set$calls, file = "genotype_data.csv", quote = F, sep = ',')
