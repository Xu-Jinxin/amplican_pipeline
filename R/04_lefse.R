#!/home/conda_envs/R-4.5.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
  library(microeco)
  library(tidyr)
  library(dplyr)
  library(tidytree)
  library(ggplot2)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript 04_lefse.R 04_lefse.conf")
} else {
  conf <- args[1]
  source(conf)
}

# source("/home/xujinxin/workflows/amplicon/R/04_lefse.conf")

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

sample.df <- read.table(file = metadata_df, sep = "\t", header = TRUE, row.names = 1)
otu.df <- read.table(file = otu_df, sep = "\t", header = TRUE, row.names = 1)
tax.df <- read.table(file = tax_df, sep = "\t", header = TRUE, row.names = 1)

otu.df <- otu.df[, sample.df$SampleName, drop = FALSE]
otu.df <- otu.df[rowSums(otu.df) > 0, , drop = FALSE]

tax.df <- tax.df[rownames(otu.df), , drop = FALSE]

tax.df <- tax.df %>%
  separate(col = Taxon, sep = ";", fill = "right",into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species) 

# filter(Kingdom != "Unassigned") %>%
# mutate(across(everything(), ~ replace_na(.x, "")))

# 创建 microtable 对象并进行 LEfSe 分析
dataset <- microtable$new(sample_table = sample.df, otu_table = otu.df, tax_table = tax.df)
lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = group, p_adjust_method = "none",
                        alpha = 0.01, lefse_subgroup = NULL)

write.table(lefse$res_diff[lefse$res_diff$LDA > lda_threshold, ], file = paste0(output, "/", lefse_table, ".tsv"), 
            sep = "\t", quote = FALSE, row.names = TRUE)

group_order <- levels(forcats::fct_inorder(sample.df[[group]]))

# 绘制 LDA 柱状图和系统发育树
p_lefse_bar <- lefse$plot_diff_bar(use_number = 1:feature_num, width = 0.8, group_order = group_order)

ggsave(paste0(output, "/", lefse_bar, ".pdf"), plot = p_lefse_bar, width = bar.width, height = bar.height)
ggsave(paste0(output, "/", lefse_bar, ".png"), plot = p_lefse_bar, width = bar.width, height = bar.height, dpi = 600)

p_lefse_cladogram <- lefse$plot_diff_cladogram(use_taxa_num = use_taxa_num, use_feature_num = feature_num, clade_label_level = clade_label_level, group_order = group_order)

ggsave(paste0(output, "/", lefse_cladogram, ".pdf"), plot = p_lefse_cladogram, width = cladogram.width, height = cladogram.height)
ggsave(paste0(output, "/", lefse_cladogram, ".png"), plot = p_lefse_cladogram, width = cladogram.width, height = cladogram.height, dpi = 600)
