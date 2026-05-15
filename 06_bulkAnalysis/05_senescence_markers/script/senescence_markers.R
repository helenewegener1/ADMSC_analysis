# Example immunomodulatory list
setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

# 3. Immunomodulatory Gene Signature Focus
# -	Provide list of known immunomodulatory genes. 
# - From DEGs, extract which immunomodulatory genes are significantly upregulated in each condition vs Control.
# -	Compare patterns of immunomodulatory gene activation between single and combined cytokine stimulations.

library(DESeq2)
library(glue)
library(stringr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(biomaRt)
library(patchwork)
library(grid)

batch_nr <- 1
# batch_nr <- 2

# ------------------------------------------------------------------------------
# Load DESeqDataSet object
# ------------------------------------------------------------------------------

dds <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))
res <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res_batch_{batch_nr}.rds'))
metadata <- dds@colData %>% as.data.frame() %>% arrange(ID) 

source("colors.R")

ann_colors <- pheatmap_colors[[glue("batch_{batch_nr}")]]

# ------------------------------------------------------------------------------
# Define markers
# ------------------------------------------------------------------------------

# Senescence markers from slides May 2026
senescence_genes_symbol <- c(
  
  # Up-regulated by genomic stress 
  "TP53", "CDKN1A", "CDKN2A", "PTEN", "TSC2",
  
  # Down-regulated by genomic stress 
  # "CDK2", "E2F1", "E2F6", "LMMB1", "MTOR",
  
  # 
  "GLB1", "IL6",
  
  # Epigenetic senescence 
  "EZH2", "DNMT1", 
  
  # SASP factors
  "COL1A1", "COL5A2", "IGFBP5", 
  
  # Actin Stress fiber senescence
  "MYLK", "CALD1", "ACTA2", "TAGLN"
  
)

# Immunomodulatory markers from slides May 2026
immunomod_genes_symbol <- c(
  
  "CD274", "IDO1", "PTGS2", "CXCL9", "CXCL10", "TGIF1", "ICAM1", "CCL2", "PTGER2"
  
)

# Prolifeation markers from slides May 2026
proliferation_genes_symbol <- c(
  
  "MKI67", "CCNB1", "CDC20", "UBE2C"
  
)

# Combine
gene_lists_symbol <- list(senescence_genes_symbol, immunomod_genes_symbol, proliferation_genes_symbol)
names(gene_lists_symbol) <- c("senescence", "immunomod", "proliferation")

# ------------------------------------------------------------------------------
# Translate gene symbols to ensembl genes
# ------------------------------------------------------------------------------

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(dds) %>% str_split_i('\\.', 1),
  mart = mart
)

# test
# grep("CD274", gene_map$hgnc_symbol, value = TRUE, ignore.case = TRUE)

# Translate gene symbols of markers to ensembl genes
gene_lists <- lapply(gene_lists_symbol, function(gl) {
  
  # gl <- senescence_genes_symbol
  genes <- gene_map %>% filter(hgnc_symbol %in% gl) %>% dplyr::select(ensembl_gene_id) %>% unlist()
  names(genes) <- gl
  return(genes)
  
})

# ==============================================================================
# Plotting
# ==============================================================================

for (markers in names(gene_lists)){

  # markers <- "senescence"
  
  # ------------------------------------------------------------------------------
  # Heatmap of gene expression across all samples.
  # ------------------------------------------------------------------------------
  
  # Get normalized expression
  vsd <- vst(dds, blind = FALSE)
  rownames(vsd) <- sub("\\..*$", "", rownames(vsd)) # Remove version numbers (everything after the dot)
  mat <- assay(vsd)[rownames(vsd) %in% (gene_lists[[markers]]), rownames(metadata)]
  
  # check order 
  table(rownames(mat) == gene_lists[[markers]])
  
  # Rename genes 
  rownames(mat) <- names(gene_lists[[markers]])
  
  # Optionally order by gene family or variance
  # mat <- mat[order(rowMeans(mat), decreasing = TRUE), ]
  
  # Plot 
  p <- pheatmap(
    mat,
    annotation_col = metadata %>% dplyr::select(ID, Donor),
    annotation_colors = ann_colors, 
    scale = "row",
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    show_colnames = TRUE, 
    angle_col = 45,
    # fontsize_row = 6, 
    fontsize_col = 8, 
    main = glue("{str_to_title(markers)} Gene Expressions"),
    silent = TRUE  # prevents drawing immediately
  )
  
  pdf(glue('06_bulkAnalysis/05_senescence_markers/plot/batch_{batch_nr}_heatmap_{markers}_all_samples.pdf'), width = 10)
  grid.newpage()
  pushViewport(viewport(x = unit(1, "cm"), width = unit(0.85, "npc"), just = "left"))
  grid.draw(p$gtable)
  popViewport()
  dev.off()
  
  # ------------------------------------------------------------------------------
  # Dot plot or boxplot of selected key genes showing per-condition expression levels
  # ------------------------------------------------------------------------------
  
  # key_genes <- list("IDO1", "NOS2", "HLA-A", "HLA-B", "HLA-C", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3")
  key_genes <- names(gene_lists[[markers]])
  
  df_plot <- mat %>% t() %>% cbind(metadata) 
  
  for (key_gene in key_genes){
    
    key_gene <- key_genes[[1]]
      
    p <- df_plot %>% 
      ggplot(aes(y = !!sym(key_gene), 
                 x = ID)) + 
      stat_summary(fun = "mean", 
                   geom = "crossbar", 
                   width = 0.4,       
                   linewidth = 0.2) +
      geom_jitter(aes(color = Donor), width = 0.1, alpha = 0.5, size = 2) + 
      scale_color_manual(values = ann_colors$Donor) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(
        title = glue('Per-condition expression level of {key_gene}'),
        subtitle = glue("{str_to_title(markers)} marker")
      ) 
    
    ggsave(glue('06_bulkAnalysis/05_senescence_markers/plot/batch_{batch_nr}_per_condition_expression_levels_{key_gene}.pdf'),
           plot = p,
           width = 10,
           height = 6)
    
  }
  
  
  # ------------------------------------------------------------------------------
  # Barplot of number of genes significantly upregulated per condition.
  # ------------------------------------------------------------------------------
  
  DEG_up_genes <- purrr::map_dfr(names(res), function(coef) {
    
    res[[coef]] %>%
      mutate(
        gene = sub("\\..*$", "", gene),  # strip version
        condition = coef
      ) %>%
      filter(gene %in% gene_lists[[markers]])
      # filter(gene %in% immunomod_genes, padj < 0.05, log2FoldChange > 1) # Only significant up-regulated hits
    
  })
  
  
  # Count number of upregulated genes per condition
  counts_per_condition <- DEG_up_genes %>%
    group_by(condition) %>%
    summarise(n_up = n())
  
  # Plot as barplot
  ggplot(counts_per_condition, 
         aes(x = condition, 
             y = n_up)) +
    geom_col() +
    theme_minimal() +
    labs(
      title = glue("Number of {markers} genes significantly upregulated per condition"),
      x = "Condition",
      y = "Number of upregulated genes"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ggsave(glue('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/batch_{batch_nr}_immunomodulatory_upregulated_per_condition.pdf'),
  #        width = 12,
  #        height = 6)
    
  
  # Table of immunomodulatory DEGs per condition.
  DEG_up_genes_list <- DEG_up_genes %>%
    group_split(condition, .keep = TRUE) %>%
    setNames(unique(DEG_up_genes$condition)) %>%
    as.list()  # <-- converts it to a plain list
  
  # Save output
  # Path to output file
  out_file <- glue("06_bulkAnalysis/05_senescence_markers/out/batch_{batch_nr}_DEG_up_{markers}.xlsx")
  
  # If file exists, delete it first (optional)
  if (file.exists(out_file)) file.remove(out_file)
  
  # Loop over your list of contrasts
  # Set Java memory BEFORE loading xlsx
  # options(java.parameters = "-Xmx4g")
  # Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home")
  library(rJava)
  library(xlsx)
  for (name in names(DEG_up_genes_list)){
    # print((name))
    write.xlsx(DEG_up_genes_list[[name]] %>% as.data.frame(), # Needs to be dataframe, cannot handle tibble
               file = out_file,
               sheetName = name,
               append = TRUE,
               row.names = FALSE)
  }

}

  