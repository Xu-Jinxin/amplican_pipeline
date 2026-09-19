#!/usr/bin/env Rscript

set.seed(123)  # 保证置换检验可重复

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
    library(vegan)
    library(dplyr)
    library(tidyverse)
    library(rdacca.hp)
    library(ggplot2)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript barplot.R barplot.conf")
} else {
  conf <- args[1]
  source(conf)
}

# setwd('/home/xujinxin/project/20260707_UK_Continetal_Shelf')
# source('/home/xujinxin/workflows/amplicon/R/05_dbRDA.conf')

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

otu <- read.table(file = otu_df, sep = "\t", header = TRUE, row.names = 1)
env <- read.table(file = env_df, sep = "\t", header = TRUE, row.names = 1)
map.df <- read.table(file = map_df, sep = "\t", header = TRUE)

otu <- as.data.frame(t(otu))
otu <- otu[map.df$SampleID, , drop = FALSE]
otu <- otu[, colSums(otu) > 0, drop = FALSE]

env <- env[map.df$SampleID, , drop = FALSE]
# 使用平均值填充给空白
num_cols <- sapply(env, is.numeric)
env[num_cols] <- lapply(env[num_cols], function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

colnames(env)
env <- env[, env_variables]

otu_hel <- decostand(otu, method = 'hellinger')
dis_bray <- vegdist(otu_hel, method = 'bray')

rda_hp <- rdacca.hp(dis_bray, env, method = 'RDA', type = 'adjR2')
rda_hp <- as.data.frame(rda_hp$Hier.part)
write.table(rda_hp, file = file.path(paste0(output, "/", "rdacca_hp_", flag, ".tsv")), sep = "\t", quote = FALSE, 
    row.names = TRUE, col.names = NA)
rda_hp

pcoa <- cmdscale(dis_bray, k = nrow(otu_hel)-1, eig = TRUE, add = TRUE)
pcoa_site <- pcoa$points

db_rda <- rda(pcoa_site~., env, scale = FALSE)
vif.cca(db_rda)

r2 <- RsquareAdj(db_rda)
rda_noadj <- r2$r.squared 
rda_adj <- r2$adj.r.squared

# exp_adj <- RsquareAdj(db_rda)$adj.r.squared * db_rda$CCA$eig/sum(db_rda$CCA$eig)
exp_adj <- RsquareAdj(db_rda)$r.squared * db_rda$CCA$eig/sum(db_rda$CCA$eig)
sum(exp_adj)
cca1_exp <- paste('RDA1:', round(exp_adj[1]*100, 2), '%')
cca2_exp <- paste('RDA2:', round(exp_adj[2]*100, 2), '%')

db_rda.scaling1 <- summary(db_rda, scaling = 1)

rda_site.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'wa')  
rda_env.scaling1 <- scores(db_rda, choices = 1:4, scaling = 1, display = 'bp')
rda_site.scaling1 <- as.data.frame(rda_site.scaling1)
rda_env.scaling1 <- as.data.frame(rda_env.scaling1)

colors = c("#1f78b4", "#33a02c", "#e31a1c")
rda_site.scaling1$Site <- map.df$Treat

ggplot(rda_site.scaling1, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Site), size = 2) +
  geom_segment(data = rda_env.scaling1, aes(x = 0, y = 0, xend = RDA1*0.9, yend = RDA2*0.9), 
               arrow = arrow(length = unit(0.1, 'cm'))) +
  geom_text(data = rda_env.scaling1, aes(RDA1, RDA2, label = rownames(rda_env.scaling1)), size = 3) +
  labs(x = cca1_exp, y = cca2_exp) +
  theme_bw()

ggsave(file.path(output, "/", output_file), width = 4.5, height = 3.5)
