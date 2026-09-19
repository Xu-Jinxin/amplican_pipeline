#!/home/conda_envs/R-4.5.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
    library(vegan)
    library(ggplot2)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript 03_beta_diversity.R 03_beta_diversity.conf")
} else {
  conf <- args[1]
  source(conf)
}

# source("/home/xujinxin/workflows/amplicon/R/03_beta_diversity.conf")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

### 自定义功能区
save_pdf <- function(file, plot = plot, width = 7, height = 7) {
  pdf(file, width = width, height = height)
  plot
  dev.off()
}

#### 执行区
taxon.df <- read.table(file = taxon_df, sep = "\t", header = TRUE, row.names = 1)
map.df <- read.table(file = map_df, sep = "\t", header = TRUE)
if (!is.null(metadata_df)) {
  metadata.df <- read.table(file = metadata_df, sep = "\t", header = TRUE)
}

#### 检查是否存在不规则样本名，按 map 表顺序进行排列
missing.id <- setdiff(map.df$SampleID, colnames(taxon.df))
if (length(missing.id) > 0) {
  stop("The following SampleIDs in map file are not found in taxon columns: ",
       paste(missing.id, collapse = ", "))
}
taxon.df <- as.data.frame(t(taxon.df[, map.df$SampleID]))
taxon.df <- taxon.df[, which(colSums(taxon.df) > 0)]
taxon.dist <- vegdist(taxon.df, method = method)

#### PERMANOVA 分析
formula <- as.formula(paste("taxon.dist ~", paste(factors, collapse = "+")))
permanova <- adonis2(formula, metadata.df, permutations = 999) 
write.table(permanova, paste0(output, "/", paste(factors, collapse = "+"), "_PERMANOVA_results.tsv"), sep = "\t", col.names = NA)

# betadisper 计算离散度
mod <- betadisper(taxon.dist, forcats::fct_inorder(map.df$Treat), type = "centroid")
upgma <- hclust(taxon.dist, method = 'average')

save_pdf(paste0(output, "/", "PERMDISP.pdf"), plot(mod), width = 6, height = 8)
save_pdf(paste0(output, "/", "PERMDISP_boxplot.pdf"), boxplot(mod), width = 10, height = 4.5)
save_pdf(paste0(output, "/", "upgma.pdf"), 
         plot(upgma, main = 'UPGMA\n(Bray-curtis distance)', sub = '', xlab = 'Sample', ylab = 'Height'), 
         width = 20, height = 4.5)

# PCoA 降维分析
pcoa <- cmdscale(taxon.dist, k = (nrow(taxon.df) - 1), eig = TRUE)
pcoa_axis <- as.data.frame(pcoa$point)
write.table(pcoa_axis, paste0(output, "/", 'pcoa_axis.tsv'), sep = '\t', col.names = NA)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

site <- pcoa_axis[ ,1:2]
site$group <- forcats::fct_inorder(map.df$Treat)

theme_custom <- theme(axis.title.x = element_text(size = x.title.size),
  axis.text.x = element_text(size = x.size, angle = x.angle, hjust = x.hjust), 
  axis.title.y = element_text(size = y.title.size), 
  axis.text.y = element_text(size = y.size, angle = y.angle, hjust = y.hjust))

p <- ggplot(data = site, aes(x = V1, y = V2, color = group)) +
    geom_point(size = 2) +
    labs(x = pcoa1, y = pcoa2) + 
    theme_bw() +
    theme_custom

ggsave(paste0(output, "/", output.file, ".pdf"), plot = p, width = width, height = height)
ggsave(paste0(output, "/", output.file, ".png"), plot = p, width = width, height = height, dpi = 600)
