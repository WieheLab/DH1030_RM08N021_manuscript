#.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.0.0')

library(stringr)
library(optparse)

merged_clonal_analysis <- read.csv("clonal_analysis/10x_merged_clones.csv")
merged_clonal_analysis_filtered <- subset(merged_clonal_analysis, Heavy_Fxnl== T & N_Heavy == 1)

heavy_rs <- read.csv("unfiltered/Heavy/VH_RecombinationSummaries.RF.txt", sep = '\t', header = T, stringsAsFactors = F)
heavy_rs_filtered <- heavy_rs[which(heavy_rs$UID %in% merged_clonal_analysis_filtered$ReadID_Heavy),]

write.table(heavy_rs_filtered, file = 'filtered/heavies.fxnl.1to1.RecombinationSummaries.txt', sep = '\t', quote = F, row.names = F)
