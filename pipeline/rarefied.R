#!/home/conda_envs/R-4.2.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
  library(vegan)
}))

asv_table <- read.table(file = '03_exported_tables/asv_table.tsv', sep = "\t", header = TRUE, row.names = 1)

set.seed(123)
rarefied_asv_table <- as.data.frame(vegan::rrarefy(t(asv_table), sample = min(colSums(asv_table))))
rarefied_asv_table <- rarefied_asv_table[, colSums(rarefied_asv_table) > 0]
rarefied_asv_table <- as.data.frame(t(rarefied_asv_table))

write.table(rarefied_asv_table, file = '03_exported_tables/rarefied_asv_table.tsv', sep = "\t", quote = FALSE, col.names = NA)
