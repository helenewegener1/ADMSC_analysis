setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

# 2. Global expression patterns
# - PCA using all expressed genes.
# - Hierarchical clustering/dendrogram.

library(tximport)
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(DESeq2)
library(ggplot2)
library(glue)
library(pheatmap)
library(RColorBrewer)
library(writexl)
# detach("package:biomaRt", unload = TRUE, character.only = TRUE)

# batch_nr <- 1
batch_nr <- 2

source("colors.R")
anno_colors <- pheatmap_colors[[glue("batch_{batch_nr}")]]

# ################################## Meta data ################################### 
# # Load and wrangle metadata
# # Batch 1: D2 & D3
# # Batch 2: D10, D11 & D12
metadata <- xlsx::read.xlsx('../data/Sample_ID.xlsx', sheetIndex = 1)

metadata_clean <- metadata %>%
  filter(!is.na(SAMPLE)) %>%
  dplyr::select(-NA.) %>%
  mutate(ID = ID %>% str_trim() %>% str_replace_all("\\s*\\+\\s*", "_"),
         ID = str_replace_all(ID, 'INFy', 'IFNy') %>% as.factor(), # InterFeroN-gamma
         ID = factor(ID, levels = c("Ctrl", "IFNy", "TNFa", "IL17a", "IFNy_IL17a",
                                    "IFNy_TNFa", "TNFa_IL17a", "IFNy_TNFa_IL17a")),
         Donor = Donor %>% str_trim() %>% as.factor(),
         Batch = case_when(Donor %in% c('2', '3') ~ '1',
                           Donor %in% c('10', '11', '12') ~ '2'
         ),
         Batch = as.factor(Batch)
  ) %>%
  filter(Batch == batch_nr) %>%
  column_to_rownames(var = 'SAMPLE')

samples <- rownames(metadata_clean)

rownames(metadata_clean) <- glue("{metadata_clean$Donor}_{metadata_clean$ID}")

# 
# # List of Kallisto output folders
# # samples <- list.files("03_kallisto_quant/out_filtered/")
# files <- file.path("03_kallisto_quant/out_filtered//", samples, "abundance.h5")
# names(files) <- samples
# 
# # Load tx2gene mapping
# tx2gene <- read.csv("05_tx2gene/out/tx2gene.tsv", sep = '\t')
# 
# # Import Kallisto results
# txi <- tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE)
# 
# # Make DESeqDataSet object
# dds <- DESeqDataSetFromTximport(txi, 
#                                 colData = metadata_clean,
#                                 design = ~ Donor + ID # Different designs tested in test_design.R
# )
# 
# colnames(dds) <- glue("{dds@colData$Donor}_{dds@colData$ID}")
# saveRDS(dds, glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))

dds <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))

################################################################################
############################### Synergy analysis ############################### 
################################################################################
# Produces terriable P-value histograms - don't think I have enough samples to make this...
# # Make DESeqDataSet object
# metadata_clean_synergy <- metadata_clean %>%
#   mutate(
#     IFNy = ifelse(str_detect(ID, 'IFNy'), 1, 0), 
#     TNFa = ifelse(str_detect(ID, 'TNFa'), 1, 0), 
#     IL17a = ifelse(str_detect(ID, 'IL17a'), 1, 0), 
#     Ctrl = ifelse(str_detect(ID, 'Ctrl'), 1, 0)
#   ) %>% 
#   select(-ID)
# 
# rownames(metadata_clean_synergy) <- samples
# 
# dds_synergy <- DESeqDataSetFromTximport(txi, 
#                                         colData = metadata_clean_synergy,
#                                         design = ~ IFNy + TNFa + IL17a + IFNy:TNFa + IFNy:IL17a + TNFa:IL17a + IFNy:TNFa:IL17a + Donor 
#                                         )
# 
# saveRDS(dds_synergy, '06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_synergy.rds')
# 
# rm(tx2gene, txi)
# 
# # Perform VST 
vsd <- vst(dds, blind = TRUE) # blind = TRUE for unsupervised PCA

################################################################################
##################################### PCA ###################################### 
################################################################################

pc_combinations <- list(
  c(1,2),
  c(1,3),
  c(2,3)
)

for (var in colnames(metadata_clean)){
  
  for (pcs in pc_combinations){
    
    # Extract PCA data
    pcaData <- plotPCA(vsd, intgroup = c(var, "Donor"), pcsToUse = c(pcs[1], pcs[2]), returnData = TRUE)
    
    # Plot with color by var and shape by Donor
    p <- ggplot(pcaData, aes(x = !!sym(glue('PC{pcs[1]}')), y = !!sym(glue('PC{pcs[2]}')), color = .data[[var]], shape = Donor)) +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = glue('PCA of {var} (PC{pcs[1]} vs PC{pcs[2]})'),
        color = glue('{var}'),
        shape = "Donor"
      ) 
    
    # Manual colors color
    if (var %in% names(anno_colors)){
      
      p <- p + scale_color_manual(values = anno_colors[[var]])
      
    }
    
    # Construct a descriptive filename
    file_name <- glue(
      '06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/batch_{batch_nr}_PCA_{var}_PC{pcs[1]}vsPC{pcs[2]}.pdf'
    )
    
    ggsave(p, filename = file_name, width = 10, height = 6)
    
  }
}

################################# PCA loadings ################################# 
# Run PCA manually using prcomp on the VST counts
# (plotPCA internally does this, but we need the rotation matrix)

# Get the VST matrix (genes x samples)
vst_mat <- assay(vsd)

# Run PCA (transposed: samples x genes, as prcomp expects observations x variables)
pca_res <- prcomp(t(vst_mat), scale. = FALSE, center = TRUE)

# Number of PCs you want to extract gene weights for
n_pcs <- 3  # adjust as needed (covers PC1, PC2, PC3 from your combinations)

# Connect to Ensembl
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract top 50 up and 50 down genes for each PC
top_genes_per_pc_list <- lapply(1:n_pcs, function(pc_i) {
  
  # Rotation/loadings for this PC (named vector, one value per gene)
  loadings <- pca_res$rotation[, pc_i]
  
  # Sort descending
  sorted_loadings <- sort(loadings, decreasing = TRUE)
  
  # Top 50 positive (up) and top 50 negative (down)
  top50_up   <- head(sorted_loadings, 100)
  top50_down <- tail(sorted_loadings, 100)
  
  df <- data.frame(
    gene     = c(names(top50_up), names(top50_down)) %>% str_split_i('\\.', 1),
    loading  = c(top50_up, top50_down),
    direction = c(rep("up", 100), rep("down", 100)),
    PC       = rep(glue("PC{pc_i}"), 200)
  )
  
  # Get mapping
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = df$gene, #%>% str_split_i('\\.', 1),
    mart = mart
  )
  
  gene_map <- gene_map %>% 
    mutate(hgnc_symbol = case_when(hgnc_symbol == "" ~ ensembl_gene_id, .default = hgnc_symbol)) %>% 
    dplyr::rename(gene = "ensembl_gene_id")
  
  df <- df %>% left_join(gene_map, by = "gene") %>% relocate(hgnc_symbol, .after = gene)
  
  return(df)
  
}) %>% setNames(glue("PC{1:3}"))

# Save to excel
write_xlsx(
  top_genes_per_pc_list,
  path = glue("06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/table/batch_{batch_nr}_top_genes_per_PC.xlsx")
)


################################################################################
###################### Hierarchical clustering/dendrogram ###################### 
################################################################################

########################## Sample distance dendogram ###########################

# vsd is your variance-stabilized DESeq2 object
vsd_mat <- assay(vsd)  # rows = genes, columns = samples

# Transpose so samples are rows
sample_dist <- dist(t(vsd_mat), method = "euclidean")

hc <- hclust(sample_dist, method = "complete")  # "complete" linkage

##################### Sample-to-sample correlation heatmap #####################

sample_cor <- cor(vsd_mat, method = "pearson")  
# sample_cor <- cor(vsd_mat, method = "spearman")

for (var in c('Donor', 'ID')){
  
  # n_colors <- metadata_clean[[var]] %>% levels() %>% length()
  
  metadata_tmp <- metadata_clean %>% arrange(!!sym(var))
  
  sample_cor <- sample_cor[rownames(metadata_tmp), rownames(metadata_tmp)]
  
  pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/batch_{batch_nr}_sample_to_sample_correlation_heatmap_{var}.pdf'),
      width = 10, 
      height = 8)
  
  print(
    pheatmap(sample_cor, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             annotation_col = metadata_tmp %>% dplyr::select(!!sym(var)),
             annotation_row = metadata_tmp %>% dplyr::select(!!sym(var)), 
             # annotation_colors = brewer.pal(8, "Set2"),
             show_rownames = TRUE,
             show_colnames = TRUE,
             scale = "row",
             main = "Sample-to-sample correlation")
  )
  
  dev.off()
  
}

pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/batch_{batch_nr}_sample_to_sample_correlation_heatmap.pdf'),
    width = 10, 
    height = 8)

print(
  pheatmap(sample_cor, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           annotation_col = metadata_clean %>% dplyr::select(ID, Donor),
           annotation_row = metadata_clean %>% dplyr::select(ID, Donor), 
           # annotation_colors = brewer.pal(8, "Set2"),
           show_rownames = TRUE,
           show_colnames = TRUE,
           scale = "row",
           main = "Sample-to-sample correlation")
)

dev.off()

##################### Sample-to-sample correlation heatmap #####################

# Optional: top 500 most variable genes since it's very heavy to run on all genes
top_var_genes <- order(rowVars(vsd_mat), decreasing = TRUE)[1:100]
sample_order <- metadata_clean %>% arrange(ID) %>% rownames()


library(biomaRt)
# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(vsd_mat)[top_var_genes] %>% str_split_i('\\.', 1),
  mart = mart
)

# If gene SYMBOL translation use that, else use ensembl ID
gene_map <- gene_map %>% mutate(hgnc_symbol = case_when(hgnc_symbol == "" ~ ensembl_gene_id, .default = hgnc_symbol))

vsd_mat_plot <- vsd_mat[top_var_genes, sample_order]

rownames(vsd_mat_plot) <- str_split_i(rownames(vsd_mat_plot), '\\.', 1)

vsd_mat_plot <- vsd_mat_plot[gene_map$ensembl_gene_id , ]

rownames(vsd_mat_plot) == gene_map$ensembl_gene_id

rownames(vsd_mat_plot) <- gene_map$hgnc_symbol

pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/batch_{batch_nr}_500_most_variable_genes_heatmap.pdf'), 
    width = 10, 
    height = 20)

pheatmap(vsd_mat_plot, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         annotation_col = metadata_clean %>% dplyr::select(ID, Donor),
         show_rownames = TRUE,
         show_colnames = TRUE,
         angle_col = 45,
         # scale = "row",
         main = "Hierarchical clustering heatmap of the 100 most variable genes") 

dev.off()






