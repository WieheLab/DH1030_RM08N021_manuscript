.libPaths('/hpc/group/dhvi/WieheLab/mb488/R_Libs/R_v4.4.3')

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-c", "--cutoff"), type="numeric", default=".05",
              help="UMI cutoff (remove chains with < cutoff*cell total umis", metavar="c"),
  make_option(c("-a", "--all_contig_annotations_csv"), type='character', default='all_contig_annotations.csv',
              help="cellranger all_contig_annotations.csv file"),
  make_option(c('-o', '--output_file'), type='character', default='all_contig_annotations_umi_filtered.csv',
              help='output file with umi filtering performed (same format as in file)'))

opt_parser = OptionParser(option_list=option_list,add_help_option=FALSE);
opt = parse_args(opt_parser);

cutoff <- opt$cutoff

annots <- read.csv(opt$all_contig_annotations_csv)
umi_sums <- annots %>% group_by(barcode, chain) %>% summarise(total = sum(umis))

annots_with_umi_sums <- merge(annots, umi_sums, by = c('barcode','chain'), all = T)
annots_with_umi_sums$cutoff_cell <- annots_with_umi_sums$total*cutoff
annots_with_umi_sums <- subset(annots_with_umi_sums, umis > cutoff_cell)

annots_with_umi_sums <- annots_with_umi_sums[,colnames(annots)]

write.csv(annots_with_umi_sums, opt$output_file, na = '', row.names = F, quote = F)
