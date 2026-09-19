#!/home/conda_envs/R-4.2.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript barplot.R barplot.conf")
} else {
  conf <- args[1]
  source(conf)
}

if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

taxon.df <- read.table(file = taxon_df, sep = "\t", header = TRUE, row.names = 1)
map.df <- read.table(file = map_df, sep = "\t", header = TRUE)

missing.id <- setdiff(map.df$SampleID, colnames(taxon.df))
if (length(missing.id) > 0) {
  stop("The following SampleIDs in map file are not found in taxon columns: ",
       paste(missing.id, collapse = ", "))
}

taxon.df <- taxon.df[, map.df$SampleID]

unclassified.raw  <- grepl("__$", rownames(taxon.df))
unclassified.df <- taxon.df[unclassified.raw, ]
classified.df <- taxon.df[!unclassified.raw, ]
rownames(classified.df) <- gsub("^[a-z]__", "", basename(gsub(";", "/", rownames(classified.df))))

classified.df$sum <- rowSums(classified.df)
classified.df <- classified.df[order(classified.df$sum, decreasing = TRUE), ]
classified.df$sum <- NULL

topN.df <- classified.df[1:topN, ]

others <- colSums(rbind(classified.df[(topN+1):dim(classified.df)[1], ], unclassified.df))
others.df <- data.frame(as.list(others))
row.names(others.df) <- "others"

plot.df <- rbind(topN.df, others.df)
taxon.order <- rownames(plot.df)

if (length(colors.list) < nrow(plot.df)) {
  stop("Insufficient number of colors, please provide at least ", nrow(plot.df), " colors.")
}

percent.df <- prop.table(as.matrix(plot.df), margin = 2) * 100
write.table(percent.df, file  = paste0(output, "/", output.file, ".txt"), col.names = NA,
            sep = "\t", quote=FALSE)

plot_long <- plot.df %>%
  tibble::rownames_to_column("Taxon") %>%
  pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Taxon = factor(Taxon, levels = taxon.order)) %>%
  mutate(Sample = factor(Sample, levels = map.df$SampleID))

# 设置自定义主题
theme_custom <- theme(axis.title.x = element_text(size = x.title.size),
  axis.text.x = element_text(size = x.size, angle = x.angle, hjust = x.hjust), 
  axis.title.y = element_text(size = y.title.size), 
  axis.text.y = element_text(size = y.size, angle = y.angle, hjust = y.hjust),
  legend.title = element_text(size = legend.title.size),
  legend.text = element_text(size = legend.text.size),
  legend.key.size = unit(legend.key.size, "inches")
)

p <- ggplot(plot_long, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_col(position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors.list) +
  labs(x = "Sample", y = "Relative abundance", fill = "Taxon") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme_custom


ggsave(paste0(output, "/", output.file, ".pdf"), plot = p, width = width, height = height)
ggsave(paste0(output, "/", output.file, ".png"), plot = p, width = width, height = height, dpi = 600)
