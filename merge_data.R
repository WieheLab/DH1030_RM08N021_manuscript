.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.4.0/')
library(Seurat)

setwd("/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/")

options(future.globals.maxSize = 1e20)

p_list <- Sys.glob('*_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
length(p_list)

for (i in p_list){
  print(i)
  load(i)
}

rm(i)
rm(p_list)

agg <- merge(control_gcLN_sub, y = c(control_gcSpleen_sub, control_memoryLN_sub, control_memoryPBMC_sub,
                                control_memorySpleen_sub, control_naiveLN_sub, control_naivePBMC_sub, control_naiveSpleen_sub, 
                                gcLN_sub, gcSpleen_sub, memoryLN_sub, memoryPBMC_sub, 
                                memorySpleen_sub, naiveLN_sub, naivePBMC_sub, naiveSpleen_sub),
             add.cell.ids = c("control_gcLN", "control_gcSpleen", "control_memoryLN", "control_memoryPBMC",
                              "control_memorySpleen", "control_naiveLN", "control_naivePBMC", "control_naiveSpleen",
                              "gcLN", "gcSpleen", "memoryLN", "memoryPBMC", 
                              "memorySpleen", "naiveLN", "naivePBMC", "naiveSpleen"),
             merge.data = T)
agg
save(agg, file = 'Merged_SeuratObject.RData')