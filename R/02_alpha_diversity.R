#!/home/conda_envs/R-4.5.3/bin/R

### 加载程序包
suppressWarnings(suppressPackageStartupMessages({
    library(vegan)
    library(forcats)
    library(ggplot2)
    library(multcompView)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error ! Usage: Rscript barplot.R barplot.conf")
} else {
  conf <- args[1]
  source(conf)
}
# source("/home/xujinxin/script/amplicon/02.alpha_diversity/alpha_index.conf")

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

Richness <- estimateR(t(taxon.df))[1, ]
Chao1 <- estimateR(t(taxon.df))[2, ]
ACE <- estimateR(t(taxon.df))[4, ]
Shannon <- diversity(t(taxon.df), index = "shannon")
Simpson <- diversity(t(taxon.df), index = "simpson")

alpha.df <- data.frame(Richness = Richness, Chao1 = Chao1, ACE = ACE, 
                       Shannon = Shannon, Simpson = Simpson)
write.table(alpha.df, file  = paste0(output, "/", output.file, ".tsv"), col.names = NA,
            sep = "\t", quote=FALSE)

alpha.df$Sample <- rownames(alpha.df)
alpha.df$Sample<- factor(alpha.df$Sample, levels = map.df$SampleID)
alpha.df$group <- forcats::fct_inorder(map.df$Treat)

theme_custom <- theme(axis.title.x = element_text(size = x.title.size),
  axis.text.x = element_text(size = x.size, angle = x.angle, hjust = x.hjust), 
  axis.title.y = element_text(size = y.title.size), 
  axis.text.y = element_text(size = y.size, angle = y.angle, hjust = y.hjust))

# 添加显著性检验
# formula <- as.formula(paste(index, "~ group"))
# fit <- aov(formula, data = alpha.df)
# tukey <- TukeyHSD(fit, "group")
# pvals <- tukey$group[, "p adj"]
# names(pvals) <- rownames(tukey$group)
# cld <- multcompLetters(pvals, threshold = 0.05)
# letter_df <- data.frame(
#   group = factor(names(cld$Letters), levels = levels(alpha.df$group)),
#   label = as.character(cld$Letters)
# )
# letter_df$y <- vapply(as.character(letter_df$group), function(g) {
#   max(alpha.df[[index]][alpha.df$group == g], na.rm = TRUE)
# }, numeric(1))

# y_range <- diff(range(alpha.df[[index]], na.rm = TRUE))
# letter_df$y <- letter_df$y + 0.05 * y_range
#  coord_cartesian(ylim = c(NA, max(letter_df$y) + 0.1 * y_range))


if ( plot_by_group ) {
  p <- ggplot(alpha.df, aes(x = group, y = .data[[index]], fill = group)) +
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot() + 
    scale_fill_manual(values = rep(color.list, 10)) +
    theme_bw() +
    theme_custom

} else {
  p <- ggplot(alpha.df, aes(x = Sample, y = .data[[index]])) +
    geom_point() +
    geom_line(aes(group = "all")) +
    theme_bw() +
    theme_custom
}

ggsave(paste0(output, "/", output.file, ".pdf"), plot = p, width = width, height = height)
ggsave(paste0(output, "/", output.file, ".png"), plot = p, width = width, height = height, dpi = 600)
