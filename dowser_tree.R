library(dowser)
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(Biostrings)
library(seqinr)
library(ggtree)

setwd('/datacommons/dhvi/10X_Genomics_data/SHIV_BG505.T33N_Infection/DH1030_tree')

# Import UCAs
bs_ucas <- readBStringSet("ucas.fasta")
df_ucas <- data.frame(name = names(bs_ucas), sequence = bs_ucas, row.names = NULL)

# Import heavies
bs_heavies <- readBStringSet("All_DH1030.heavies.masked.fasta")
df_heavies <- data.frame(name = names(bs_heavies), sequence = bs_heavies, row.names = NULL)
df_heavies <- df_heavies |> 
  mutate(
    chain = "heavy"
  )

# Import kappas
bs_kappas <- readBStringSet("All_DH1030.kappas.fasta")
df_kappas <- data.frame(name = names(bs_kappas), sequence =bs_kappas, row.names = NULL)
df_kappas <- df_kappas |> 
  mutate(
    chain = "kappa"
  )

df_sequences <- bind_rows(df_heavies, df_kappas) #|> 
# filter(!grepl("I", name))

#### Create Air Compartible Dataframe

sample_id <- "DH1030"
clone_id <- 1

df_airr <- data.frame(sequence_id = df_sequences$name, sequence = df_sequences$sequence, rev_comp = NA, productive = NA, v_call = NA, d_call = NA, j_call = NA, sequence_alignment = df_sequences$sequence, germline_alignment = NA, junction = NA, junction_aa = NA, v_cigar = NA, d_cigar = NA, j_cigar = NA, stop_codon = NA, vj_in_frame = NA, locus = NA , c_call = NA, junction_length = NA, np1_length = NA, np2_length = NA, v_sequence_start = NA, v_sequence_end = NA, v_germline_start = NA, v_germline_end = NA, d_sequence_start = NA, d_sequence_end = NA, d_germline_start = NA, d_germline_end = NA, j_sequence_start = NA,  j_sequence_end = NA, j_germline_start = NA, j_germline_end = NA, v_score = NA, v_identity = NA , v_support = NA, d_score = NA, d_identity = NA, d_support = NA, j_score = NA, j_identity = NA, j_support = NA, fwr1 = NA, fwr2 = NA, fwr3 = NA, fwr4 = NA, cdr1 = NA, cdr2 = NA, cdr3=NA, clone_id = clone_id, germline_v_call = NA, germline_d_call = NA, germline_j_call = NA, cell_id = df_sequences$name, sample_id = sample_id, chain = df_sequences$chain)

#### Refine AIRR df
## Add additional fields necessary
light_chain = "IGK"

df_airr <- df_airr %>% mutate(
  sequence_id = if_else(!chain == "heavy", paste0(sequence_id, "_", 2), sequence_id),
  seq_len = nchar(sequence),
  germline_alignment = if_else(chain == "heavy", df_ucas[1,2], df_ucas[2,2]),
  locus = if_else(chain == "heavy", "IGH", light_chain)
)

## Take a quick sneak peak at the contents of the file.
glimpse(df_airr)

#### Resolve light chains 
df_airr_res <- resolveLightChains(df_airr)
print(df_airr_res$clone_subgroup)

VH_RecombinationSummaries_clones <- read.csv("cloanalyst/Heavy/All_DH1030.heavies.masked.RecombinationSummaries.txt", sep = '\t')

VH_calls <- data.frame(VH_RecombinationSummaries_clones) %>% 
  select(UID, VGene, DGene, JGene, CDR3Length) |> 
  dplyr::rename(sequence_id  = UID, v_call = VGene, d_call = DGene, j_call = JGene, cdr3_length = CDR3Length) |> 
  mutate(junction_length = cdr3_length+6)

VK_RecombinationSummaries_clones <- read.csv('cloanalyst/Kappa/All_DH1030.kappas.RecombinationSummaries.txt', sep = '\t')

VK_calls <- VK_RecombinationSummaries_clones |> 
  select(UID, VGene, DGene, JGene, CDR3Length) |> 
  dplyr::rename(sequence_id  = UID, v_call = VGene, d_call = DGene, j_call = JGene, cdr3_length = CDR3Length) |> 
  mutate(
    sequence_id = paste0(sequence_id, "_2"),
    junction_length = cdr3_length+6
  )

RS_Orig <- bind_rows(VH_calls, VK_calls)

df_airr_res[ ,c("v_call", "d_call", "j_call","cdr3_length", "productive",  "junction_length")] <- list(NULL)

## Now read the columns with data
df_airr_res <-  df_airr_res |> 
  left_join(RS_Orig, by = "sequence_id")

# read in DH numbers and add new ones

dh_nums <- read.csv('/datacommons/dhvi/mb488/dowser/DH1030_tree/DH1030_nums_Pos_56_Res.csv')
dh_nums$Sequence_ID <- gsub('VDJ_RM08N021wk105', 'VDJ_RM08N021wk105_G1', dh_nums$Sequence_ID)
dh_nums$Sequence_ID <- gsub('LIBRAseqrun3_RM08N021wk105', 'LIBRAseqRun3_RM08N021wk105_G1', dh_nums$Sequence_ID)
dh_nums$Sequence_ID <- gsub('LIBRAseqrun2_RM08N021wk105', 'LIBRAseqRun2_RM08N021wk105_G1', dh_nums$Sequence_ID)

all_mems <- VH_calls %>% select(sequence_id)
all_mems$sequence_id2 <- gsub('ShawHIVRAD', 'PhenotypicBEAM', all_mems$sequence_id)

all_mems <- merge(all_mems, dh_nums, by.x = 'sequence_id2', by.y = 'Sequence_ID', all = T)

all_mems2 <- subset(all_mems, !is.na(sequence_id))
all_mems3 <- all_mems2 %>% subset(is.na(DH_Num)) %>% mutate(DH_Num = paste0('DH1030.', 200:205))
all_mems4 <- subset(all_mems2, !is.na(DH_Num))

all_mems_final <- smartbind(all_mems3, all_mems4)

write.csv(all_mems_final, file = 'All_DH1030_mems_DHnums_P56.csv', row.names = F, quote = F)

# get new members residue 
seqs <- read.fasta('All_DH1030.heavies.aa.aligned.fasta', forceDNAtolower = F, as.string = T)
seqs_df <- data.frame(sequence_id = names(seqs), Seq = as.character(seqs))

seqs_df <- seqs_df %>%
  mutate(Seq = str_split(Seq, "")) %>%  # Split into list of characters
  unnest_wider(Seq, names_sep = "") %>% data.frame()
seqs_df_sub <- select(seqs_df, sequence_id, Seq56)

new_members <- all_mems_final %>% subset(is.na(Position_56))
seqs_df_new <- subset(seqs_df_sub, sequence_id %in% new_members$sequence_id)
colnames(seqs_df_new) <- c('sequence_id', 'Position_56')

all_mems_final <- merge(all_mems_final, seqs_df_new, by = 'sequence_id', all = T)
all_mems_final$Position_56_out <- ifelse(is.na(all_mems_final$Position_56.x), all_mems_final$Position_56.y, all_mems_final$Position_56.x)
all_mems_final <- select(all_mems_final, sequence_id, DH_Num, Position_56_out)

#### Format clones At this point, each row is sequence and each column is some information about that sequence. However, for tree inference, we need to reformat our data such that, each row is a clone and each column is information about that clone. The function will collapse identical sequences together with the same trait, if the traits differ e.g. sampling date, those sequences will remain distinct. Numerical fields will be added together when collapsing occurs. You can specify a minimum number of sequences that a clone should have to be retained, otherwise it is dropped. The final output is a *AirClone* where each row is a clone, with clone_id, number of sequences as well as other additional information.
## Sort out lengths

clones <- formatClones(df_airr_res, chain="HL", nproc=1, collapse = FALSE, use_regions = F,
                       split_light = T, minseq = 3, germ = "germline_alignment")
print(clones)

### Build the tree

# Building maximum likelihood trees with multiple partitions using IgPhyML 
# Only the newest version of IgPhyML supports this option
# exec here is set to IgPhyML position in the Docker image.
## The omega warning can be avoided if you set the omega variable to "e,e" when calling getTrees (when making paired trees) and that should go away. With the sequences that are shared, it is likely the "N"s are added to make the sequence length a multiple of three. 
gt_raxml <- getTrees(clones, build="raxml", nproc=1, rm_temp = F, dir = "temp",
                     exec="/datacommons/dhvi/jes183/software/raxml-ng/build/bin/raxml-ng", collapse = F)
  
gt_raxml$data[[1]]@data <- merge(gt_raxml$data[[1]]@data, all_mems_final, by = 'sequence_id', all = T)
gt_raxml$data[[1]]@data$Position_56_out2 <- ifelse(gt_raxml$data[[1]]@data$Position_56_out %in% c('G', 'R'), gt_raxml$data[[1]]@data$Position_56_out, 'Other' )

tree <- plotTrees(gt_raxml, tips = 'Position_56_out2', guide_title = 'Position 56 Residue')[[1]] + 
  geom_tiplab(offset =  0, size = 3, mapping = aes(label = DH_Num)) + 
  ggtitle("DH1030 Lineage Tree (Dowser)") +
  coord_cartesian(clip="off") +  ggplot2::xlim(0, .1) +
  scale_color_manual(values = c('green','red', 'blue'), limits = c('G', 'Other', 'R'), breaks = c('G', 'R', 'Other')) 

pdf("DH1030_tree_new.raxml.DH_Nums.Position56.pdf", height = 30, width = 8)
tree
dev.off()

png("DH1030_tree_new.raxml.DH_Nums.Position56.png", height = 30, width = 8, units = 'in', res = 600)
tree
dev.off()

# get BEAM scores immunogenetics

beam <- read.csv('../Reverse_BEAM/CR7.2/RM08N021/wk105/Combined_Clonal_Analysis/')



all_seqs <- getAllSeqs(gt_raxml)
all_heavy_seqs <- subset(all_seqs, locus == 'IGH')
all_kappa_seqs <- subset(all_seqs, locus == 'IGK')

write.fasta(sequences = as.list(all_heavy_seqs$sequence), names = all_heavy_seqs$node_id, 
            file.out = 'all_heavies_with_intermediates_germline.raxml.fasta', nbchar = 10000, as.string = T)

write.fasta(sequences = as.list(all_kappa_seqs$sequence), names = all_kappa_seqs$node_id, 
            file.out = 'all_kappas_with_intermediates_germline.raxml.fasta', nbchar = 10000, as.string = T)

# disambiguate intermediates
p_gt_raxml <- plotTrees(gt_raxml, node_nums = T, labelsize = 3)[[1]]

p_gt_raxml_data <- p_gt_raxml$data
intermediate_node_labs <- p_gt_raxml_data |> 
  filter(isTip == FALSE) |> 
  select(node, label)
## Always set the node number here:

uca_node_num <- min(intermediate_node_labs$node)

## IGH UCA:
df_raxml_UCA_IGH = data.frame(node = uca_node_num, sequence = getNodeSeq(gt_raxml, node=uca_node_num, clone="1_1")[[1]])
df_raxml_UCA_IGH <- df_raxml_UCA_IGH |> 
  mutate(
    sequence_clean = gsub("[.]", "", str_sub(sequence, 1, nchar(sequence)-2)),
    seq_len = nchar(sequence_clean)
  ) |> 
  select(seq_len, node, sequence_clean) |> 
  dplyr::rename(sequence = sequence_clean)

## IGH intermediates: Ancestral probabilities of intermediates inferred by dowser
dowser_uca_ints_AP <- read.delim("temp/sample_1_1_1_asr.raxml.ancestralProbs")

dowser_uca_ints_AP_hc <- dowser_uca_ints_AP 

dowser_uca_ints_AP_hc_updated <- dowser_uca_ints_AP |> 
  rowwise() |> 
  mutate(
    site_disambig = if_else(Site %in% c("A","C","T","G"), State,
                            sapply(str_split(c("p_A", "p_C", "p_G", "p_T")[which.max(across(p_A:p_T))], "_", 2),`[`, 2)
    )
  ) |> 
  ungroup()

dowser_uca_ints_AP_hc_updated_seq <- dowser_uca_ints_AP_hc_updated |> 
  select(Node, site_disambig) |> 
  group_by(Node) |> 
  summarise(sequence = gsub(",","", toString(site_disambig), fixed = T)) |> 
  ungroup() |> 
  mutate(sequence = gsub("[[:blank:]]", "", sequence),
         sequence = str_sub(sequence,1, nchar(sequence)-2),
         seq_len = nchar(sequence))
dowser_uca_ints_AP_hc_updated_seq <- dowser_uca_ints_AP_hc_updated_seq |> 
  inner_join(intermediate_node_labs, by = c("Node"="label")) |> 
  filter(!node == uca_node_num) 

## Add UCA to the top
df_raxml_UCA_intermediates_IGH <- bind_rows(df_raxml_UCA_IGH ,dowser_uca_ints_AP_hc_updated_seq)

write.fasta(sequences = as.list(df_raxml_UCA_intermediates_IGH$sequence),
            names = if_else(df_raxml_UCA_intermediates_IGH$node == uca_node_num, 
                            paste0("i", df_raxml_UCA_intermediates_IGH$node, "_UCA"), 
                            paste0("i", df_raxml_UCA_intermediates_IGH$node)),
            file.out = "pool/dowser/DH570_raxml_intermediates_dowser_IGH.fasta")

# get intermediates
all_sequences = getAllSeqs(gt_raxml)
intermediates <- subset(all_sequences, grepl(pattern = 'Node', x = all_sequences$node_id) & locus == 'IGH')
intermediates_kappa <- subset(all_sequences, grepl(pattern = 'Node', x = all_sequences$node_id) & locus == 'IGK')

write.fasta(sequences = as.list(intermediates$sequence), names = intermediates$node_id, file.out = 'intermediates.fasta', nbchar = 100000, as.string = T)
write.fasta(sequences = as.list(intermediates_kappa$sequence), names = intermediates_kappa$node_id, file.out = 'intermediates.kappa.fasta', nbchar = 100000, as.string = T)

png("DH1030_tree_new.raxml.DH_Nums.NodeNums.png", height = 50, width = 12, units = 'in', res = 600)
plotTrees(gt_raxml, node_nums = T)[[1]] + 
  geom_tiplab(offset =  0, size = 3, mapping = aes(label = DH_Num)) + 
  ggtitle("DH1030 Lineage Tree (Dowser)") +
  coord_cartesian(clip="off") +  ggplot2::xlim(0, .1) 
dev.off()


# intermediates

## Node names
p_gt_raxml_data <- tree$data
intermediate_node_labs <- p_gt_raxml_data |> 
  filter(isTip == FALSE) |> 
  select(node, label)

## IGH intermediates: Ancestral probabilities of intermediates inferred by dowser
dowser_uca_ints_AP <- read.delim("temp/sample_1_1_1_asr.raxml.ancestralProbs")

dowser_uca_ints_AP_hc <- dowser_uca_ints_AP 
dowser_uca_ints_AP_hc_updated <- dowser_uca_ints_AP_hc |> 
  rowwise() |> 
  mutate(
    site_disambig = if_else(Site %in% c("A","C","T","G"), State,
                            sapply(str_split(c("p_A", "p_C", "p_G", "p_T")[which.max(across(p_A:p_T))], "_", 2),`[`, 2)
    )
  ) |> 
  ungroup()

dowser_uca_ints_AP_hc_updated_seq <- dowser_uca_ints_AP_hc_updated |> 
  select(Node, site_disambig) |> 
  group_by(Node) |> 
  summarise(sequence = gsub(",","", toString(site_disambig), fixed = T)) |> 
  ungroup() |> 
  mutate(sequence = gsub("[[:blank:]]", "", sequence),
         sequence = str_sub(sequence,1, nchar(sequence)-2),
         seq_len = nchar(sequence))
dowser_uca_ints_AP_hc_updated_seq <- dowser_uca_ints_AP_hc_updated_seq |> 
  inner_join(intermediate_node_labs, by = c("Node"="label")) |> 
  filter(!node == uca_node_num) 

## Add UCA to the top
df_raxml_UCA_intermediates_IGH <- bind_rows(df_raxml_UCA_IGH ,dowser_uca_ints_AP_hc_updated_seq)

write.fasta(sequences = as.list(df_raxml_UCA_intermediates_IGH$sequence),
            names = if_else(df_raxml_UCA_intermediates_IGH$node == uca_node_num, 
                            paste0("i", df_raxml_UCA_intermediates_IGH$node, "_UCA"), 
                            paste0("i", df_raxml_UCA_intermediates_IGH$node)),
            file.out = "intermediates.IGH.fasta")

## IGL UCA:
df_raxml_UCA_IGL = data.frame(node = uca_node_num, sequence = getNodeSeq(gt_raxml_CN, node=uca_node_num, clone="1_1")[[2]])
df_raxml_UCA_IGL <- df_raxml_UCA_IGL |> 
  mutate(
    sequence_clean = gsub("[.]", "", str_sub(sequence, 1, nchar(sequence)-2)),
    Node = "Node5",
    seq_len = nchar(sequence_clean)
  ) |> 
  select(Node, seq_len, node, sequence_clean) |> 
  dplyr::rename(sequence = sequence_clean)

dowser_uca_ints_AP_lc <- dowser_uca_ints_AP |> 
  filter(Part == 2)
dowser_uca_ints_AP_lc_updated <- dowser_uca_ints_AP_lc |> 
  rowwise() |> 
  mutate(
    site_disambig = if_else(Site %in% c("A","C","T","G"), State,
                            sapply(str_split(c("p_A", "p_C", "p_G", "p_T")[which.max(across(p_A:p_T))], "_", 2),`[`, 2)
    )
  ) |> 
  ungroup()

dowser_uca_ints_AP_lc_updated_seq <- dowser_uca_ints_AP_lc_updated |> 
  select(Node, site_disambig) |> 
  group_by(Node) |> 
  summarise(sequence = gsub(",","", toString(site_disambig), fixed = T)) |> 
  ungroup() |> 
  mutate(sequence = gsub("[[:blank:]]", "", sequence),
         sequence = str_sub(sequence,1, nchar(sequence)-2),
         seq_len = nchar(sequence))
dowser_uca_ints_AP_lc_updated_seq <- dowser_uca_ints_AP_lc_updated_seq |> 
  inner_join(intermediate_node_labs, by = c("Node"="label")) |> 
  filter(!node == uca_node_num) 

## Add UCA to the top
df_raxml_UCA_intermediates_IGL <- bind_rows(df_raxml_UCA_IGL ,dowser_uca_ints_AP_lc_updated_seq)

write.fasta(sequences = as.list(df_raxml_UCA_intermediates_IGL$sequence),
            names = if_else(df_raxml_UCA_intermediates_IGL$node == uca_node_num, 
                            paste0("i", df_raxml_UCA_intermediates_IGL$node, "_UCA"), 
                            paste0("i", df_raxml_UCA_intermediates_IGL$node)),
            file.out = "dowser/DH270_raxml_intermediates_dowser_IGL.fasta")
