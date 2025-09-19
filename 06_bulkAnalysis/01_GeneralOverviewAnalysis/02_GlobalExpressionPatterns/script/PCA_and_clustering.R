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
detach("package:biomaRt", unload = TRUE, character.only = TRUE)


# List of Kallisto output folders
samples <- list.files("03_kallisto_quant/out")
files <- file.path("03_kallisto_quant/out/", samples, "abundance.h5")
names(files) <- samples

# Load tx2gene mapping
tx2gene <- read.csv("05_tx2gene/out/tx2gene.tsv", sep = '\t')

# Import Kallisto results
txi <- tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE)

# Load and wrangle metadata
# Batch 1: D2 & D3
# Batch 2: D10, D11 & D12
metadata <- xlsx::read.xlsx('../data/Sample_ID.xlsx', sheetIndex = 1)

################################################################################
######################## cytokine condition vs Control #########################
################################################################################

metadata_clean <- metadata %>%
  filter(SAMPLE %in% samples) %>%
  dplyr::select(-NA.) %>%
  mutate(ID = ID %>% str_trim() %>% str_replace_all("\\s*\\+\\s*", "_"),
         ID = str_replace_all(ID, 'INFy', 'IFNy') %>% as.factor(), # InterFeroN-gamma
         Donor = Donor %>% str_trim() %>% as.factor(),
         Batch = case_when(Donor %in% c('2', '3') ~ '1',
                           Donor %in% c('10', '11', '12') ~ '2'
                           ),
         Batch = as.factor(Batch)
         ) %>%
  column_to_rownames(var = 'SAMPLE')

rownames(metadata_clean) <- samples

# Make DESeqDataSet object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = metadata_clean,
                                design = ~ Donor + ID # Different designs tested in test_design.R
                                )

saveRDS(dds, '06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds.rds')

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
    
    # Construct a descriptive filename
    file_name <- glue(
      '06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/PCA_{var}_PC{pcs[1]}vsPC{pcs[2]}.pdf'
    )
    
    ggsave(p, filename = file_name, width = 10, height = 6)
  }
}

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
  
  pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/sample_to_sample_correlation_heatmap_{var}.pdf'),
      width = 10, 
      height = 8)
  
  print(
    pheatmap(sample_cor, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             annotation_col = metadata_tmp %>% select(!!sym(var)),
             annotation_row = metadata_tmp %>% select(!!sym(var)), 
             # annotation_colors = brewer.pal(8, "Set2"),
             show_rownames = TRUE,
             show_colnames = TRUE,
             scale = "row",
             main = "Sample-to-sample correlation")
  )
  
  dev.off()
  
}

pdf(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/sample_to_sample_correlation_heatmap.pdf'),
    width = 10, 
    height = 8)

print(
  pheatmap(sample_cor, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           annotation_col = metadata_clean %>% select(ID, Donor),
           annotation_row = metadata_clean %>% select(ID, Donor), 
           # annotation_colors = brewer.pal(8, "Set2"),
           show_rownames = TRUE,
           show_colnames = TRUE,
           scale = "row",
           main = "Sample-to-sample correlation")
)

dev.off()

##################### Sample-to-sample correlation heatmap #####################

# Optional: top 500 most variable genes since it's very heavy to run on all genes
top_var_genes <- order(rowVars(vsd_mat), decreasing = TRUE)[1:500]

pdf('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/plot/500_most_variable_genes_heatmap.pdf', 
    width = 10, 
    height = 10)

pheatmap(vsd_mat[top_var_genes, ], 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = metadata_clean[, c("ID", "Donor"), drop=FALSE],
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "row",
         main = "Hierarchical clustering heatmap of the 500 most variable genes") 

dev.off()






