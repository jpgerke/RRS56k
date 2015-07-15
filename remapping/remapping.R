rm(list=ls())
require(dplyr)
require(stringr)
load("../data/new_full_set.RData")

fromchris <- read.delim("../dataforupload/55kv2_sorted.tab", stringsAsFactors = F, header = T)

remap <- read.delim("../dataforupload/cluster_b_remapped_calls.txt", stringsAsFactors = F, header = T)

df <- data.frame(Chrom = new_full_set$LG,
                 Position = new_full_set$position,
                 Name = as.character(new_full_set$markernames),
                 Code = as.character(new_full_set$markercodes),
                 Gen = new_full_set$F2map,
                 stringsAsFactors = F)

rm(new_full_set)

original = select(fromchris, Name, Code=dbSNP, Chrom_Original=chr.v2, Position_Original=pos.v2) 

test = left_join(df, original, by=c("Name", "Code")) %>% 
  mutate(Chrom_Original = as.integer(str_replace(Chrom_Original, "chr", ""))) %>%
  mutate(Position_Original = as.integer(Position_Original))

#how many markers have moved chromosomes?
sum(test$Chrom != test$Chrom_Original)
#[1] 15
samechrom = filter(test, Chrom == Chrom_Original)
sum(abs(samechrom$Position - samechrom$Position_Original>0) )
#[1] 1585

moved = c()
for (mychrom in 1:10) {
  currchrom = filter(samechrom, Chrom==mychrom)
  tots = sum(abs(currchrom$Position - currchrom$Position_Original>0) )
  moved = c(moved, tots)
}

moved
#[1]   20   41  221    2  106   90   39    2 1063    1


#look at chr9 out of curiosity
chrom9 = filter(samechrom, Chrom==9) 
chrom9_gen = left_join(chrom9, df)
moved_9 = filter(chrom9_gen, Position != Position_Original) %>%
  mutate(difference = abs(Position - Position_Original))
ggplot(data=moved_9) + aes(x=Position, y=Gen) + geom_point() #+ geom_hline(yintercept=138.4)
ggplot(data=moved_9) + aes(x=Position, y=Position_Original) + geom_point()
