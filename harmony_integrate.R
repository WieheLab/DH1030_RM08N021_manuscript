.libPaths('/datacommons/dhvi/mb488/R_Libs/R_v4.4.0/')
library(Seurat)

setwd("/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/2_integration/harmony_integration/")

load(file = '../../1_preprocessing/Merged_SeuratObject.RData')

# usually remove Ig genes, however only 15 so many are probably unannotated
# all_genes <- data.frame(genes = rownames(agg))
# ig_genes <- all_genes[grepl('IG[HKL]V|IG[KL]J|IG[KL]C', all_genes$genes),]
# all_non_ig_genes <- all_genes[!(all_genes$genes %in% ig_genes),]
# 
# agg
# agg <- subset(agg, features = all_non_ig_genes)

agg <- NormalizeData(agg)
agg <- FindVariableFeatures(agg)
#all_genes <- rownames(agg)
#agg <- ScaleData(agg, features = all_genes)
agg <- ScaleData(agg)
agg <- RunPCA(agg, npcs = 30, verbose = F)

pdf('ElbowPlot.pdf', height = 6, width = 8)
ElbowPlot(agg, ndims = 30)
dev.off()

agg <- IntegrateLayers(object = agg, method = HarmonyIntegration, 
                       orig.reduction = "pca", 
                       new.reduction = "harmony", verbose = TRUE)

agg[["RNA"]] <- JoinLayers(agg[["RNA"]])

#agg <- FindNeighbors(agg, reduction = "integrated.cca", dims = 1:30)
agg <- RunUMAP(agg, dims = 1:30, reduction = "harmony")

save(agg, file = 'Agg.harmony_integrated.RData')