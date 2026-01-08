library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
memorySpleen_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_spleen/GEX/RM08N021_MemSpleen_Agg/outs/count/filtered_feature_bc_matrix/')
memorySpleen <- CreateSeuratObject(counts = memorySpleen_data, project = "memorySpleen")

# basic qc
#memorySpleen[["percent.mt"]] <- PercentageFeatureSet(memorySpleen, pattern = "^MT-")
VlnPlot(memorySpleen, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

memorySpleen <- subset(memorySpleen, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
memorySpleen <- NormalizeData(memorySpleen)
memorySpleen <- FindVariableFeatures(memorySpleen, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memorySpleen)
memorySpleen <- ScaleData(memorySpleen, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/memory_spleen/BEAM/ShawHIVRAD_RM08N021wk105MemSpleen_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
vdj_new <- vdj
vdj_new$barcode2 <- gsub('-[1-8]', '', vdj_new$barcode)
vdj_new$lane <- gsub('[AGCT]{16}-', '', vdj_new$barcode)

vdj_new$NewLane <- case_when(vdj_new$lane == '1' ~ '2',
                             vdj_new$lane == '2' ~ '1',
                             vdj_new$lane == '3' ~ '3',
                             vdj_new$lane == '4' ~ '4')

vdj_new$NewBarcode <- paste(vdj_new$barcode2, as.character(vdj_new$NewLane), sep = '-')
memorySpleen@meta.data$VDJ_data <- ifelse(rownames(memorySpleen@meta.data) %in% vdj_new$NewBarcode, T, F)
table(memorySpleen@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(memorySpleen), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(memorySpleen), ref = ref, labels = ref$label.fine)
memorySpleen[["SingleR.labels.main"]] <- singler.labels.main$labels
memorySpleen[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj_new) <- vdj_new$NewBarcode
vdj_new <- dplyr::select(vdj_new, -NewBarcode, -NewLane, -barcode, -barcode2, -lane)
memorySpleen <- AddMetaData(memorySpleen, vdj_new)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
memorySpleen_sub <- subset(memorySpleen, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
memorySpleen_sub <- NormalizeData(memorySpleen_sub)
memorySpleen_sub <- FindVariableFeatures(memorySpleen_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(memorySpleen_sub)
memorySpleen_sub <- ScaleData(memorySpleen_sub, features = all.genes)

save(memorySpleen_sub, file = 'memorySpleen_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(memorySpleen_sub@meta.data)

# memorySpleen_sub <- subset(memorySpleen, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


