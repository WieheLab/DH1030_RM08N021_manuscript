library(Seurat)
library(SingleR)
library(SeuratObject)
library(ggplot2)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/1_preprocessing/')

# read in data and create seurat object
naivePBMC_data <- Read10X('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_PBMC/GEX/cellranger/RM08N021_NaivePBMC_Agg/outs/count/filtered_feature_bc_matrix/')
naivePBMC <- CreateSeuratObject(counts = naivePBMC_data, project = "naivePBMC")

# basic qc
#naivePBMC[["percent.mt"]] <- PercentageFeatureSet(naivePBMC, pattern = "^MT-")
VlnPlot(naivePBMC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = .1)

naivePBMC <- subset(naivePBMC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# normalize, scale data
naivePBMC <- NormalizeData(naivePBMC)
naivePBMC <- FindVariableFeatures(naivePBMC, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naivePBMC)
naivePBMC <- ScaleData(naivePBMC, features = all.genes)

# read in fxnl 1:1 vdj data with BEAM data
vdj <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/naive_PBMC/BEAM/ShawHIVRAD_RM08N021wk105NaivePBMC_merged_clonal_analysis.fxnl_1to1.BEAM_Scores.UMIs.Seqs.csv')
vdj <- subset(vdj, !(is.na(vdj$BG505_Mut11b_SOSIP_1_BEAM_Score)))
vdj_new <- vdj
vdj_new$barcode2 <- gsub('-[1-8]', '', vdj_new$barcode)
vdj_new$lane <- gsub('[AGCT]{16}-', '', vdj_new$barcode)

vdj_new$NewLane <- case_when(vdj_new$lane == '1' ~ '3',
                             vdj_new$lane == '2' ~ '1',
                             vdj_new$lane == '3' ~ '2')

vdj_new$NewBarcode <- paste(vdj_new$barcode2, as.character(vdj_new$NewLane), sep = '-')
naivePBMC@meta.data$VDJ_data <- ifelse(rownames(naivePBMC@meta.data) %in% vdj_new$NewBarcode, T, F)
table(naivePBMC@meta.data$VDJ_data)

# run singleR to be able to subset only B cells
ref <- celldex::MonacoImmuneData()
bcell_subsets <- c('Exhausted B cells', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Plasmablasts')
singler.labels.main <- SingleR(as.SingleCellExperiment(naivePBMC), ref = ref, labels = ref$label.main)
singler.labels.fine <- SingleR(as.SingleCellExperiment(naivePBMC), ref = ref, labels = ref$label.fine)
naivePBMC[["SingleR.labels.main"]] <- singler.labels.main$labels
naivePBMC[["SingleR.labels.fine"]] <- singler.labels.fine$labels

# add vdj metadata
rownames(vdj_new) <- vdj_new$NewBarcode
vdj_new <- dplyr::select(vdj_new, -NewBarcode, -NewLane, -barcode, -barcode2, -lane)
naivePBMC <- AddMetaData(naivePBMC, vdj_new)

# subset for only B cells (main and fine labels) with fxnl 1:1 vdj data
naivePBMC_sub <- subset(naivePBMC, SingleR.labels.main == 'B cells' & VDJ_data == T & SingleR.labels.fine %in% bcell_subsets)

# normalize, scale for integration with all samples
naivePBMC_sub <- NormalizeData(naivePBMC_sub)
naivePBMC_sub <- FindVariableFeatures(naivePBMC_sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(naivePBMC_sub)
naivePBMC_sub <- ScaleData(naivePBMC_sub, features = all.genes)

save(naivePBMC_sub, file = 'naivePBMC_Bcells_wFxnl1to1VDJandBEAM_SeuratObject.RData')
nrow(naivePBMC_sub@meta.data)

# naivePBMC_sub <- subset(naivePBMC, seurat_clusters != 3 & SingleR.labels.main == 'B cells')


