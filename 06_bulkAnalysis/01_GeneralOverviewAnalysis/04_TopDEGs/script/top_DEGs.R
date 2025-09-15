setwd('~/Documents/projects/project_ADMSC/analysis/')

# 4.	Top DEGs
# o	Identify top 50 upregulated and top 50 downregulated genes for each comparison.
# •	Heatmap of top 50 DEGs for selected comparisons.
# •	Volcano plots for each comparison.
# •	Barplots of number of DEGs per comparison (up vs down).

library(dplyr)
library(pheatmap)
library(glue)
library(DESeq2)
library(stringr)
library(ggplot2)

# Load data
dds <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds.rds')
res <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res.rds')
# resLFC <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/resLFC.rds')

# filter top 50 up and down regulated genes 
res_filtered <- list()
resLFC_filtered <- list()

for (name in names(res)){
  
  res_tmp <- res[[name]]
  
  # Sort by log2FoldChange (descending) for top upregulated genes
  top_up <- head(res_tmp[order(res_tmp$log2FoldChange, decreasing = TRUE), ], 50)
  
  # Sort by log2FoldChange (ascending) for top downregulated genes
  top_down <- head(res_tmp[order(res_tmp$log2FoldChange, decreasing = FALSE), ], 50)
  
  # Optionally, combine them for export
  top_combined <- rbind(top_up, top_down)
  
  res_filtered[[name]] <- top_combined
  
  # Clean up 
  rm(top_up, top_down, top_combined)
  
  # # Log fold change shrinkage
  # resLFC_tmp <- resLFC[[name]]
  # top_up <- head(resLFC_tmp[order(resLFC_tmp$log2FoldChange, decreasing = TRUE), ], 50)
  # top_down <- head(resLFC_tmp[order(resLFC_tmp$log2FoldChange, decreasing = FALSE), ], 50)
  # top_combined <- rbind(top_up, top_down)
  # 
  # resLFC_filtered[[name]] <- top_combined
  # 
  # # Clean up 
  # rm(top_up, top_down, top_combined)
  
}

# Save output
# Path to output file
out_file <- "06_bulkAnalysis/01_GeneralOverviewAnalysis/04_TopDEGs/out/res_filtered.xlsx"

# If file exists, delete it first (optional)
if (file.exists(out_file)) file.remove(out_file)

# Loop over your list of contrasts
# Set Java memory BEFORE loading xlsx
# options(java.parameters = "-Xmx4g")
# Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home")
library(rJava)
library(xlsx)
for (name in names(res_filtered)){
  write.xlsx(res_filtered[[name]],
             file = out_file,
             sheetName = name,
             append = TRUE,
             row.names = FALSE)
}

# Path to output file
out_file <- "06_bulkAnalysis/01_GeneralOverviewAnalysis/04_TopDEGs/out/resLFC_filtered.xlsx"

# If file exists, delete it first (optional)
if (file.exists(out_file)) file.remove(out_file)

# Loop over your list of contrasts
# Set Java memory BEFORE loading xlsx
# options(java.parameters = "-Xmx4g")
# Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home")
library(rJava)
library(xlsx)
for (name in names(res_filtered)){
  write.xlsx(resLFC_filtered[[name]],
             file = out_file,
             sheetName = name,
             append = TRUE,
             row.names = FALSE)
}

################################################################################
#################################### Plots ##################################### 
################################################################################

meta_data <- dds@colData

for (comparison in names(res)){
  
  # comparison <- names(res)[1]
  
  res_tmp <- res[[comparison]]
  
  # 🦊	Heatmap of top 50 DEGs for selected comparisons.
  
  # Top 50 DEGs (based on padj)
  res_arranged <- res_tmp[order(res_tmp$padj),]
  top_genes <- head(res_arranged, n = 50) %>% dplyr::select(gene) %>% unlist()
  
  # Extract variance-stabilized counts
  vsd <- vst(dds, blind = TRUE)
  vsd_mat <- assay(vsd)
  
  # Subset vsd for top genes
  mat_top <- vsd_mat[top_genes, ]
  
  # Assess samples of conditions 
  IDs <- str_split(comparison, '_vs_') %>% unlist()
  samples_of_condition <- meta_data %>% as.data.frame() %>% filter(ID %in% IDs) 
  
  # Subset mat_top to samples of condition 
  mat_top_subset <- mat_top[, colnames(mat_top) %in% rownames(samples_of_condition)] 
  mat_top_subset <- mat_top_subset[, order(samples_of_condition$ID)]
  
  pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/04_TopDEGs/plot/top_50_DEGs_heatmap_{comparison}.pdf'))
  
  p <- pheatmap(mat_top_subset,
                annotation_col = samples_of_condition %>% select(ID),
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                fontsize_row = 6,
                main = glue("Top 50 DEGs in {comparison}"))
  
  print(p)
  
  dev.off()
  
  # 🦊	Volcano plots for each comparison.
  df_volcano <- res_tmp %>%
    mutate(
      neglog10padj = -log10(padj),
      significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                           TRUE, FALSE)
    )
  
  df_volcano %>% 
    ggplot(aes(x = log2FoldChange, 
               y = neglog10padj)) +
    geom_point(aes(color = significant), alpha = 0.7) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_bw() +
    labs(
      title = glue("Volcano Plot: {comparison}"),
      subtitle = "DEGs: padj < 0.05 & abs(log2FoldChange) > 1", 
      x = "log2 Fold Change",
      y = "-log10(padj)",
      color = "DEGs"
    )
  
  ggsave(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/04_TopDEGs/plot/volcano_{comparison}.pdf'), 
         width = 10, 
         height = 7)
  
}

# 🦊	Barplots of number of DEGs per comparison (up vs down).
deg_summary <- lapply(names(res), function(comparison) {
  res_tmp <- res[[comparison]]
  res_tmp <- res_tmp[!is.na(res_tmp$padj), ]
  
  up <- sum(res_tmp$padj < 0.05 & res_tmp$log2FoldChange > 0)
  down <- sum(res_tmp$padj < 0.05 & res_tmp$log2FoldChange < 0)
  
  data.frame(
    comparison = comparison,
    Direction = c("Up", "Down"),
    count = c(up, down)
  )
}) %>% bind_rows()

deg_summary %>% 
  ggplot(aes(x = comparison, 
             y = count, 
             fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c('grey20', 'grey80')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of DEGs per Comparison",
    subtitle = "Down: padj < 0.05 & log2FoldChange < 0 \nUp: padj < 0.05 & log2FoldChange > 0 ",
    x = "Comparison",
    y = "Number of DEGs"
  ) + 
  theme(
    plot.subtitle = element_text(size = 8)  
  )

ggsave('06_bulkAnalysis/01_GeneralOverviewAnalysis/04_TopDEGs/plot/N_DEGs_all_comparisons.pdf', 
       width = 10, 
       height = 7)















