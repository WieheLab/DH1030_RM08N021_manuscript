files <- Sys.glob('*/outs/metrics_summary.csv')

out <- data.frame()

for (i in files){
  m <- read.csv(i)
  cells <- as.numeric(gsub(",", "",  m$Estimated.Number.of.Cells))
  mrpc <- as.numeric(gsub(",", "",  m$Mean.Read.Pairs.per.Cell))
  reads <- mrpc*cells
  print(paste0( i, ': ', cells, ' Cells, ', mrpc, ' Mean Reads Per Cell'))
  out <- rbind(out, data.frame(Cells = cells, TotalReads = reads))
}

print(paste0('Total Cells: ', sum(out$Cells)))
print(paste0('Mean Reads Per Cell: ', round(sum(out$TotalReads)/sum(out$Cells), digits = 1)))

