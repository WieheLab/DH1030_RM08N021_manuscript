library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
gcSpleen_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/gc_spleen/GEX/RM08N021_GCSpleen_Agg/outs/count/filtered_feature_bc_matrix/')
gcSpleen <- CreateSeuratObject(counts = gcSpleen_data, project = "gcSpleen")

# basic qc
#gcSpleen[["percent.mt"]] <- PercentageFeatureSet(gcSpleen, pattern = "^MT-")
VlnPlot(gcSpleen, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

gcSpleen <- subset(gcSpleen, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
gcSpleen <- NormalizeData(gcSpleen)
gcSpleen <- FindVariableFeatures(gcSpleen, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gcSpleen)
gcSpleen <- ScaleData(gcSpleen, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/gc_spleen/BEAM/ShawHIVRAD_RM08N021wk105GCSpleen_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
gcSpleen@meta.data$VDJ_data <- ifelse(rownames(gcSpleen@meta.data) %in% vdj$barcode, T, F)
table(gcSpleen@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(gcSpleen), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(gcSpleen), ref = ref, labels = ref$label.fine)
gcSpleen[["SingleR.labels.main"]] <- singler.labels.main$labels
gcSpleen[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj) <- vdj$barcode
vdj <- dplyr::select(vdj, -barcode)
gcSpleen <- AddMetaData(gcSpleen, vdj)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
gcSpleen_sub <- subset(gcSpleen, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
gcSpleen_sub <- NormalizeData(gcSpleen_sub)
gcSpleen_sub <- FindVariableFeatures(gcSpleen_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gcSpleen_sub)
gcSpleen_sub <- ScaleData(gcSpleen_sub, features = all.genes)

save(gcSpleen_sub, file = 'gcSpleen_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(gcSpleen_sub@meta.data)

# gcSpleen_sub <- subset(gcSpleen, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


