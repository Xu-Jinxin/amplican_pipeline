#!/usr/bin/env.df Rscript

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
    library(vegan)
    library(dplyr)
    library(tidyverse)
    library(rdacca.hp)
    library(ggplot2)
}))

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("Error ! Usage: Rscript barplot.R barplot.conf")
# } else {
#   conf <- args[1]
#   source(conf)
# }
set.seed(123)  # 保证置换检验可重复

setwd('/home/xujinxin/project/20260707_UK_Continetal_Shelf')
source('/home/xujinxin/workflows/amplicon/R/06_mantel.conf')

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

otu.df <- read.table(file = otu_df, sep = "\t", header = TRUE, row.names = 1)
env.df <- read.table(file = env_df, sep = "\t", header = TRUE, row.names = 1)
map.df <- read.table(file = map_df, sep = "\t", header = TRUE)

otu.df <- as.data.frame(t(otu.df))
otu.df <- otu.df[map.df$SampleID, , drop = FALSE]
otu.df <- otu.df[, colSums(otu.df) > 0, drop = FALSE]

env.df <- env.df[map.df$SampleID, , drop = FALSE]
# 使用平均值填充给空白
num_cols <- sapply(env.df, is.numeric)
env.df[num_cols] <- lapply(env.df[num_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

dist.otu <- vegdist(otu.df, method = "bray")
dist.otu.list <- as.vector(dist.otu)

mantel.result <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(mantel.result) <- c("Factor", "Mantel_r", "p_value")

for (factor in colnames(env.df)) {
    dist.env <- vegdist(env.df[[factor]], method = "euclidean")
    otu_env <- mantel(dist.otu, dist.env, method = 'spearman', permutations = 999, na.rm = TRUE)
    dist.env.list <- as.vector(dist.env)
    ggplot(data = data.frame(dist.otu.list, dist.env.list), aes(x = dist.env.list, y = dist.otu.list)) +
        geom_point(size = 2, alpha = 0.8) +
        geom_smooth(method = "lm", formula = 'y ~ x', se = FALSE, color = "blue") +
        labs(title = paste("Mantel Test: ", factor),
             x = paste(factor, "Distance"),
             y = "Microbial Distance") +
        annotate("text", x = max(dist.env.list) * 0.7, y = max(dist.otu.list) * 0.9,
                 label = paste("Mantel r:", round(otu_env$statistic, 3), "\n",
                               "p-value:", otu_env$signif),
                 hjust = 0)

    ggsave(filename = paste0(output, "/", factor, "_mantel_plot.png"), width = 6, height = 4, dpi = 300)
    ggsave(filename = paste0(output, "/", factor, "_mantel_plot.pdf"), width = 6, height = 4)
    
    
    mantel.result.temp <- data.frame(factor = factor, mantel_r = otu_env$statistic, p_value = otu_env$signif)
    mantel.result <- rbind(mantel.result, mantel.result.temp)
}

write.table(mantel.result, file = paste0(output, "/mantel_results.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


