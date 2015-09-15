rm(list=ls())
library(ggplot2)
library(cowplot)


load("plot_BSCB1.rds")
load("plot_BSSS.rds")

combined = plot_grid(plot_BSSS, plot_BSCB1, labels=c("A", "B"))
save_plot("../fig_S1_combined.pdf", combined,
          ncol = 2, 
          nrow = 1, 
          base_aspect_ratio = 1.1
)

combined = plot_grid(plot_BSSS, plot_BSCB1, labels=c("A", "B"))
save_plot("../fig_5.jpg", combined,
          ncol = 2, 
          nrow = 1, 
          base_aspect_ratio = 1.1
)
