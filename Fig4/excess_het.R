
rm(list=ls())
library(ggplot2)
library(dplyr)

simresults = read.table("../data/2cm_sim_collated_results.txt", header=T, stringsAsFactors = F) %>%
  mutate(Genlength = genstop - genstart, Physlength = (phystop-phystart)/1000000)

excess_BSSS <- simresults %>% filter(BSSS_upperP < 0.001)
excess_BSCB <- simresults %>% filter(BSCB_upperP < 0.001)

