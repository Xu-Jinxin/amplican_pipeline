#!/home/conda_envs/R-4.2.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(picante)
    library(iCAMP)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript barplot.R barplot.conf")
} else {
  conf <- args[1]
  source(conf)
}
# source('/home/xujinxin/workflows/amplicon/R/icamp.conf')

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}


map.df <- read.table(file = map_df, header = TRUE)
tree <- read.tree(file = tree_file)
taxon.df <- read.table(file = taxon_df, sep = '\t', header = TRUE, row.names = 1)
taxon.df  <- as.data.frame(t(taxon.df)) # 读取 OTU 表，以样本作为行名，OTU作为列名

factors <- unique(map.df[[group]])

for (factor in factors) {
    taxon.df2 <- taxon.df [which(rownames(taxon.df) %in% map.df[which(map.df[[group]] %in% factor), ]$SampleID), ]
    taxon.df2 <- taxon.df2[, which(colSums(taxon.df2) > 0)]
    prune_tree <- prune.sample(taxon.df2, tree) 
    pd <- cophenetic(prune_tree)
    qpen.out <- qpen(comm = taxon.df2, pd = pd, sig.bNTI = 2, sig.rc = 0.95, 
                 rand.time = 999, nworker = 8)
    ratio <- qpen.out$ratio
    ratio <- melt(ratio, variable.name = 'process', value.name = 'prop')
    ratio <- subset(ratio, process != 'num.pair')
    ratio$process <- sapply(ratio$process, function(x){gsub( '\\.', ' ', x)})
    write.table(ratio, file = paste0(output, '/', factor, '_iCAMP.txt'), col.names = NA,
            sep = "\t", quote=FALSE)
}