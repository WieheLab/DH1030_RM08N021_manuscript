library(Seurat)
library(ggplot2)
library(dplyr)

setwd("/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/3_visualization")

load(file = '/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/Analysis/GEX/Seurat/2_integration/harmony_integration/Agg.harmony_integrated.RData')

agg@meta.data$SampleType <- gsub('control_', '', agg@meta.data$orig.ident)
agg@meta.data$Monkey <- ifelse(grepl('control', x = agg@meta.data$orig.ident), 'Control', 'RM08N021')
agg@meta.data$Compartment <- case_when(grepl(pattern = 'gc', agg@meta.data$SampleType) ~ 'Germinal Center',
                                       grepl(pattern = 'memory', agg@meta.data$SampleType) ~ 'Memory',
                                       grepl(pattern = 'naive', agg@meta.data$SampleType) ~ 'Naive')
agg@meta.data$Source <- case_when(grepl(pattern = 'LN', agg@meta.data$SampleType) ~ 'Lymph Node',
                                       grepl(pattern = 'PBMC', agg@meta.data$SampleType) ~ 'PBMC',
                                       grepl(pattern = 'Spleen', agg@meta.data$SampleType) ~ 'Spleen')

png('Orig.Ident.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'orig.ident', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Sample') 
dev.off()

png('Orig.Ident.DimPlot.Split.png', height = 18, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'orig.ident', split.by = 'orig.ident', raster = F, shuffle = T, ncol = 4) + 
  labs(title = 'Reverse BEAM\nUMAP by Sample')  + NoLegend()
dev.off()

png('Monkey.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey') 
dev.off()

png('Monkey.DimPlot.Split.png', height = 9.5, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', raster = F, shuffle = T, split.by = 'Monkey') + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey')  + NoLegend()
dev.off()

set.seed(4293)
min_cells <- min(table(agg@meta.data$Monkey))
agg@meta.data$rownames <- rownames(agg@meta.data)
cells_sample <- agg@meta.data %>% 
  group_by(Monkey) %>%
  slice_sample(n = min_cells) %>% 
  data.frame() %>% select(rownames)

png('Monkey_Sampled.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', raster = F, shuffle = T, cells = cells_sample$rownames) + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey (Subsampled)') 
dev.off()

png('Monkey_Sampled.DimPlot.Split.png', height = 9.5, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', raster = F, shuffle = T, split.by = 'Monkey', cells = cells_sample$rownames) + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey (Subsampled)')  + NoLegend()
dev.off()

png('SampleType.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'SampleType', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Sample Type') 
dev.off()

png('SampleType.DimPlot.Split.png', height = 9, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'SampleType', split.by = 'SampleType', raster = F, shuffle = T, ncol = 4) + 
  labs(title = 'Reverse BEAM\nUMAP by Sample Type')  + NoLegend()
dev.off()

png('Monkey_SampleType.DimPlot.Split.png', height = 9, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', split.by = 'SampleType', raster = F, shuffle = T, ncol = 4) + 
  labs(title = 'Reverse BEAM\nUMAP by Sample Type & Monkey') 
dev.off()

png('SampleType_Monkey.DimPlot.Split.png', height = 9.5, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'SampleType', raster = F, shuffle = T, split.by = 'Monkey') + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey & Sample Type')  
dev.off()

png('Compartment.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'Compartment', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Compartment') 
dev.off()

png('Compartment.DimPlot.Split.png', height = 8, width = 20, units = 'in', res = 600)
DimPlot(agg, group.by = 'Compartment', split.by = 'Compartment', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by Compartment')  + NoLegend()
dev.off()

png('Monkey_Compartment.DimPlot.Split.png', height = 7, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', split.by = 'Compartment', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by Compartment & Monkey') 
dev.off()

png('Compartment_Monkey.DimPlot.Split.png', height = 9.5, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'Compartment', raster = F, shuffle = T, split.by = 'Monkey') + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey & Compartment')  
dev.off()

png('Source.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'Source', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Source') 
dev.off()

png('Source.DimPlot.Split.png', height = 8, width = 20, units = 'in', res = 600)
DimPlot(agg, group.by = 'Source', split.by = 'Source', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by Source')  + NoLegend()
dev.off()

png('Monkey_Source.DimPlot.Split.png', height = 7, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'Monkey', split.by = 'Source', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by Source & Monkey') 
dev.off()

png('Source_Monkey.DimPlot.Split.png', height = 9.5, width = 18.5, units = 'in', res = 600)
DimPlot(agg, group.by = 'Source', raster = F, shuffle = T, split.by = 'Monkey') + 
  labs(title = 'Reverse BEAM\nUMAP by Monkey & Source')  
dev.off()

png('SingleR.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'SingleR.labels.fine', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by SingleR Classification') 
dev.off()

png('SingleR.DimPlot.Split.png', height = 12, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'SingleR.labels.fine', split.by = 'SingleR.labels.fine', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by SingleR Classification') + NoLegend()
dev.off()

png('Source_SingleR.DimPlot.Split.png', height = 8, width = 19, units = 'in', res = 600)
DimPlot(agg, group.by = 'SingleR.labels.fine', split.by = 'Source', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by SingleR Classification & Source')
dev.off()

png('Compartment_SingleR.DimPlot.Split.png', height = 8, width = 19, units = 'in', res = 600)
DimPlot(agg, group.by = 'SingleR.labels.fine', split.by = 'Compartment', raster = F, shuffle = T, ncol = 3) + 
  labs(title = 'Reverse BEAM\nUMAP by SingleR Classification & Compartment') 
dev.off()

png('SampleType_SingleR.DimPlot.Split.png', height = 12, width = 21, units = 'in', res = 600)
DimPlot(agg, group.by = 'SingleR.labels.fine', split.by = 'SampleType', raster = F, shuffle = T, ncol = 4) + 
  labs(title = 'Reverse BEAM\nUMAP by SingleR Classification & Sample Type') 
dev.off()

features <- c('MS4A1', 'TFRC', 'CD38', )

colnames(agg@meta.data)

# cell cmopartment markers

Idents(agg) <- agg@meta.data$Compartment
#compartment_markers <- FindAllMarkers(agg, assay = 'RNA', )
# write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.csv', quote = F, row.names = F)

compartment_markers <- read.csv(file = 'RevBEAM_CompartmentMarkers.csv')

avg_expression <- AverageExpression(agg, 
                                    assays = 'RNA', 
                                    features = unique(compartment_markers$gene), 
                                    group.by = 'Compartment', 
                                    verbose = T)

colnames(avg_expression$RNA) <- paste(gsub(' ', '', colnames(avg_expression$RNA)), '_AverageExpression', sep = '')
avg_expression_df <- data.frame(avg_expression$RNA)
avg_expression_df$gene <- row.names(avg_expression$RNA)

compartment_markers <- merge(compartment_markers, avg_expression_df, by = 'gene')

compartment_markers <- select(compartment_markers, gene, cluster, avg_log2FC, pct.1, pct.2, 
                              GerminalCenter_AverageExpression:Naive_AverageExpression, 
                              p_val, p_val_adj)

compartment_markers <- compartment_markers[order(compartment_markers$cluster, abs(compartment_markers$avg_log2FC), decreasing = T),]
write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.with_AvgExpression.csv', quote = F, row.names = F)

# compartment markers RM08 only

rm08 <- subset(agg, Monkey == 'RM08N021')
Idents(rm08) <- rm08@meta.data$Compartment
compartment_markers <- FindAllMarkers(rm08, assay = 'RNA', logfc.threshold = 0, )
write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.RM08N021_only.csv', quote = F, row.names = F)

compartment_markers <- read.csv(file = 'RevBEAM_CompartmentMarkers.RM08N021_only.csv')

avg_expression <- AverageExpression(rm08, 
                                    assays = 'RNA', 
                                    features = unique(compartment_markers$gene), 
                                    group.by = 'Compartment', 
                                    verbose = T)

colnames(avg_expression$RNA) <- paste(gsub(' ', '', colnames(avg_expression$RNA)), '_AverageExpression', sep = '')
avg_expression_df <- data.frame(avg_expression$RNA)
avg_expression_df$gene <- row.names(avg_expression$RNA)

compartment_markers <- merge(compartment_markers, avg_expression_df, by = 'gene')

compartment_markers <- select(compartment_markers, gene, cluster, avg_log2FC, pct.1, pct.2, 
                              GerminalCenter_AverageExpression:Naive_AverageExpression, 
                              p_val, p_val_adj)

compartment_markers <- compartment_markers[order(compartment_markers$cluster, abs(compartment_markers$avg_log2FC), decreasing = T),]
write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.with_AvgExpression.RM08N021_only.csv', quote = F, row.names = F)

# compartment markers naive RM only

naiveOnly <- subset(agg, Monkey == 'Control')
Idents(naiveOnly) <- naiveOnly@meta.data$Compartment
compartment_markers <- FindAllMarkers(naiveOnly, assay = 'RNA')
write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.naive_only.csv', quote = F, row.names = F)

compartment_markers <- read.csv(file = 'RevBEAM_CompartmentMarkers.naive_only.csv')

avg_expression <- AverageExpression(naiveOnly, 
                                    assays = 'RNA', 
                                    features = unique(compartment_markers$gene), 
                                    group.by = 'Compartment', 
                                    verbose = T)

colnames(avg_expression$RNA) <- paste(gsub(' ', '', colnames(avg_expression$RNA)), '_AverageExpression', sep = '')
avg_expression_df <- data.frame(avg_expression$RNA)
avg_expression_df$gene <- row.names(avg_expression$RNA)

compartment_markers <- merge(compartment_markers, avg_expression_df, by = 'gene')

compartment_markers <- select(compartment_markers, gene, cluster, avg_log2FC, pct.1, pct.2, 
                              GerminalCenter_AverageExpression:Naive_AverageExpression, 
                              p_val, p_val_adj)

compartment_markers <- compartment_markers[order(compartment_markers$cluster, abs(compartment_markers$avg_log2FC), decreasing = T),]
write.csv(compartment_markers, file = 'RevBEAM_CompartmentMarkers.with_AvgExpression.naive_only.csv', quote = F, row.names = F)

# compartment DEG marker heatmap
avg_expression <- AverageExpression(agg, 
                                    assays = 'RNA', 
                                    features = unique(compartment_markers$gene), 
                                    group.by = 'Compartment', 
                                    verbose = T, return.seurat = T)

top10_genes <- compartment_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = abs(avg_log2FC), n = 10, with_ties = FALSE) %>%
  ungroup()

png('Avg_Expression_Top10_Compartment.Heatmap.png', height = 6, width = 6, units = 'in', res = 600)
DoHeatmap(avg_expression, features = top10_genes$gene, slot = 'data', draw.lines = F)
dev.off()

# cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
agg <- CellCycleScoring(agg, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

png('CellCyclePhase.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'Phase', raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nUMAP by Cell Cycle Phase') 
dev.off()

cell_cycle_compartment_freqs <- agg@meta.data %>% 
  group_by(Compartment, Phase) %>%
  summarize(n = n()) %>%
  mutate(Frequency = n / sum(n)) %>%
  data.frame()

png('CellCyclePhase.BarPlot.png', height = 6, width = 6, units = 'in', res = 600)
ggplot(cell_cycle_compartment_freqs, aes(x = Compartment, y = Frequency, fill = Phase)) + 
  geom_bar(stat = 'identity', color = 'black') + 
  labs(title = 'Reverse BEAM\nCell Cycle Phase Frequency by Compartment')  + 
  theme(plot.title = element_text(hjust = .5))
dev.off()


# clustering 

agg <- FindNeighbors(agg, reduction = "harmony", dims = 1:30)
for (i in seq(.05, 1.1, by = .05)){
  agg <- FindClusters(agg, resolution = i)
  print(i)
}

x <- clustree(agg, prefix = "RNA_snn_res.")

pdf('clustree.pdf', height = 10, width = 12)
x
dev.off()

pdf('MYC_clustree.pdf', height = 8, width = 10)
clustree(agg, node_colour = "MYC", node_colour_aggr = "mean")
dev.off()

sc3_stability <- x$data %>% group_by(RNA_snn_res.) %>% summarize(mean = mean(sc3_stability))

pdf('sc3_stability.pdf', height = 3, width = 6)
ggplot(sc3_stability, aes(x = as.numeric(as.character(RNA_snn_res.)), y = mean)) + 
  geom_point() + 
  geom_line(group = 1) + 
  labs(title = 'SC3 Stability', x = 'Resolution', y = 'Mean SC3 Stability') + 
  geom_hline(yintercept = max(sc3_stability$mean), color = 'red', linetype = 'dashed') + 
  theme(plot.title = element_text(hjust = .5)) + 
  scale_x_continuous(breaks = seq(0, 1.1, by = .1))
dev.off()

agg <- FindClusters(agg, resolution = .5)

png('Cluster.DimPlot.png', height = 10, width = 12, units = 'in', res = 600)
DimPlot(agg, group.by = 'seurat_clusters', label = T, raster = F)  + 
  labs(title = 'Reverse BEAM\nUMAP by Cluster') 
dev.off()

png('Cluster.DimPlot.Split.png',  height = 18, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'seurat_clusters', raster = F, split.by = 'seurat_clusters', ncol = 4)  + 
  labs(title = 'Reverse BEAM\nUMAP by Cluster') + NoLegend()
dev.off()

png('Cluster_Compartment.DimPlot.Split.png',  height = 8, width = 21, units = 'in', res = 600)
DimPlot(agg, group.by = 'seurat_clusters', raster = F, split.by = 'Compartment', ncol = 3, label = T)  + 
  labs(title = 'Reverse BEAM\nUMAP by Cluster and Compartment') 
dev.off()

png('Compartment_Cluster.DimPlot.Split.png',  height = 18, width = 18, units = 'in', res = 600)
DimPlot(agg, group.by = 'Compartment', raster = F, split.by = 'seurat_clusters', ncol = 4, label = F)  + 
  labs(title = 'Reverse BEAM\nUMAP by Compartment and Cluster') 
dev.off()

cluster_compartment_freqs <- agg@meta.data %>% 
  group_by(seurat_clusters, Compartment) %>%
  summarize(n = n()) %>%
  mutate(Frequency = n / sum(n)) %>%
  data.frame()

png('Cluster_Compartment.BarPlot.png', height = 6, width = 10, units = 'in', res = 600)
ggplot(cluster_compartment_freqs, aes(x = seurat_clusters, y = Frequency, fill = Compartment)) + 
  geom_bar(stat = 'identity', color = 'black') + 
  labs(title = 'Reverse BEAM\nCompartment Frequency by Cluster')  + 
  theme(plot.title = element_text(hjust = .5))
dev.off()

# cluster DEGs - run in tmux
Idents(agg) <- agg@meta.data$RNA_snn_res.0.5
cluster_markers <- FindAllMarkers(agg, assay = 'RNA')

cluster_markers <- select(cluster_markers, gene, cluster, avg_log2FC, pct.1, pct.2, 
                              p_val, p_val_adj)

write.csv(cluster_markers, file = 'RevBEAM_ClusterMarkers.csv', quote = F, row.names = F)


igd1 <- FetchData(agg, vars = c('ENSMMUG00000061441', 'ENSMMUG00000057110'))
igd1$sum <- (igd1$ENSMMUG00000061441 > 0) + (igd1$ENSMMUG00000057110>0)

rm08 <- ScaleData(rm08, features = rownames(rm08))

features <- c('ENSMMUG00000061441', 'ENSMMUG00000057110', 'IGHM', 'FOXP1', 'FCER2', 'CR2')
pdf('Naive_BCell_Genes.DotPlot.pdf', height = 4, width = 8)
DotPlot(rm08, features = features, group.by = 'Compartment', scale = T) + 
  labs(title = 'Conventional Naive B Cell Genes') + 
  theme(plot.title = element_text(hjust = .5), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


features <- c('ENSMMUG00000061441', 'ENSMMUG00000057110', 'IGHM', 'FOXP1', 'FCER2', 'CR2')
pdf('Naive_BCell_Genes.not_scaled.DotPlot.pdf', height = 4, width = 8)
DotPlot(rm08, features = features, group.by = 'Compartment', scale = F) + 
  labs(title = 'Conventional Naive B Cell Genes') + 
  theme(plot.title = element_text(hjust = .5), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

features <- c('ENSMMUG00000061441', 'ENSMMUG00000057110', 'IGHM', 'FOXP1', 'FCER2', 'CR2')
rm08 <- ScaleData(rm08, features = features)
avg_scale_expression <- AverageExpression(rm08, layer = 'scale.data', features = features, group.by = 'Compartment')
write.csv(data.frame(avg_scale_expression$RNA), file = 'Naive_B_Scaled_AvgExpression.csv', row.names = T)

features <- c('BCL6', 'MKI67', 'AICDA', 'S1PR2', 'JCHAIN')
rm08 <- ScaleData(rm08, features = features)
avg_scale_expression <- AverageExpression(rm08, layer = 'scale.data', features = features, group.by = 'Compartment')
write.csv(data.frame(avg_scale_expression$RNA), file = 'GC_B_Scaled_AvgExpression.csv', row.names = T)

features <- c('ITGAX', 'TBX21', 'CD86', 'TNFRSF13B', 'ENSMMUG00000015202', 'ENSMMUG00000040771')
rm08 <- ScaleData(rm08, features = features)
avg_scale_expression <- AverageExpression(rm08, layer = 'scale.data', features = features, group.by = 'Compartment')
write.csv(data.frame(avg_scale_expression$RNA), file = 'Memory_B_Scaled_AvgExpression.csv', row.names = T)

load('RM08.SeuratObject.RData')

rm08_ln <- subset(rm08, Source == 'Lymph Node')
rm(rm08)

rm08_ln[["RNA"]] <- split(rm08_ln[["RNA"]], f = rm08_ln$orig.ident)
rm08_ln <- NormalizeData(rm08_ln)
rm08_ln <- FindVariableFeatures(rm08_ln)
rm08_ln <- ScaleData(rm08_ln)
rm08_ln <- RunPCA(rm08_ln)
ElbowPlot(rm08_ln, ndims = 30)
rm08_ln <- IntegrateLayers(object = rm08_ln, method = HarmonyIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "harmony", verbose = TRUE)

# re-join layers after integration
rm08_ln[["RNA"]] <- JoinLayers(rm08_ln[["RNA"]])
rm08_ln <- RunUMAP(rm08_ln, dims = 1:24, reduction = "harmony")

png('RM08N021_LNonly.Compartment.DimPlot.png', height = 8, width = 11, units = 'in', res = 600)
DimPlot(rm08_ln, raster = F, shuffle = T) + 
  labs(title = 'Reverse BEAM\nRM08N021 Lymph Node\nUMAP by Compartment') +
  theme(plot.title = element_text(hjust = .5))
dev.off()

dh1030_mems <- read.csv('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/Reverse_BEAM/CR7.2/RM08N021/wk105/DH1030/RevBEAM_DH1030_new_members.BEAM.with_new_members.csv')
dh1030_mems$BCell_Compartment <- gsub('GCLN', 'gcLN', dh1030_mems$BCell_Compartment)
dh1030_mems$BCell_Compartment <- gsub('Mem', 'memory', dh1030_mems$BCell_Compartment)
dh1030_gcln <- subset(dh1030_mems, BCell_Compartment %in% c('gcLN'))
dh1030_gcln$lane <- gsub('[AGCT]{16}-', '', dh1030_gcln$barcode)
dh1030_gcln$lane_fixed <- case_when(dh1030_gcln$lane == '3' ~'5',
                                    dh1030_gcln$lane == '4' ~'3',
                                    dh1030_gcln$lane == '5' ~'6',
                                    dh1030_gcln$lane == '6' ~'4',
                                    dh1030_gcln$lane == '1' ~'1',
                                    dh1030_gcln$lane == '2' ~'2',
                                    dh1030_gcln$lane == '7' ~'7')
dh1030_gcln$barcode2 <- gsub('[1-8]', '', dh1030_gcln$barcode)
dh1030_gcln$new_id <- paste(dh1030_gcln$BCell_Compartment, paste0(dh1030_gcln$barcode2, dh1030_gcln$lane_fixed), sep = '_')
dh1030_memln <- subset(dh1030_mems, BCell_Compartment %in% c('memoryLN'))

dh1030_ln <- smartbind(dh1030_gcln, dh1030_memln)
dh1030_ln$new_id2 <- ifelse(is.na(dh1030_ln$new_id), paste(dh1030_ln$BCell_Compartment, dh1030_ln$barcode, sep = '_'), dh1030_ln$new_id)

dh1030_cells <- dh1030_ln$new_id2
table(dh1030_cells %in% rownames(rm08_ln@meta.data))
dh1030_cells[!(dh1030_cells %in% rownames(rm08_ln@meta.data))]

png('RM08N021_LNonly.DH1030_Members.DimPlot.png', height = 8, width = 11, units = 'in', res = 600)
DimPlot(rm08_ln, cells.highlight = dh1030_cells, raster = F, sizes.highlight = 3) + 
  scale_color_manual(labels = c("Unselected", "DH1030 Clonal Members"), values = c("grey", "red")) + 
  labs(title = 'Reverse BEAM\nRM08N021 Lymph Node\nUMAP with DH1030 Clonal Members') + 
  theme(plot.title = element_text(hjust = .5))
dev.off()
  

dh1030_mems_sub$Position_56_Group <- ifelse(dh1030_mems_sub$Position_56_out %in% c('G', 'R'), dh1030_mems_sub$Position_56_out, 'Other')
