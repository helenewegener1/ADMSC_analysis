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

# Load DESeqDataSet object
dds <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds.rds')
res <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res.rds')
metadata <- dds@colData %>% as.data.frame() %>% arrange(ID) 

# Pulled out the ass - find real ones... 
# immunomod_genes_symbol <- c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "TNFRSF9", "CD80", "CD86", "IL10", "IDO1")
immunomod_genes_symbol <- c("CXCL10", "IL6")

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(dds) %>% str_split_i('\\.', 1),
  mart = mart
)

immunomod_genes <- gene_map %>% filter(hgnc_symbol %in% immunomod_genes_symbol) %>% dplyr::select(ensembl_gene_id) %>% unlist()
names(immunomod_genes) <- immunomod_genes_symbol

# 🐸	Heatmap of immunomodulatory gene expression across all samples.

# Get normalized expression
vsd <- vst(dds)
rownames(vsd) <- sub("\\..*$", "", rownames(vsd)) # Remove version numbers (everything after the dot)
mat <- assay(vsd)[rownames(vsd) %in% immunomod_genes, rownames(metadata)]

# Optionally order by gene family or variance
# mat <- mat[order(rowMeans(mat), decreasing = TRUE), ]

pdf('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/heatmap_immunomodulatory_all_samples.pdf')
pheatmap(mat,
         annotation_col = metadata %>% dplyr::select(ID),
         scale = "row",
         show_rownames = TRUE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         fontsize_row = 6,
         main = "Immunomodulatory Gene Expression Across Samples")
dev.off()

# 🐸	Dot plot or boxplot of selected key genes showing per-condition expression levels.
key_gene <- immunomod_genes[2]

df_plot <- mat %>% t() %>% cbind(metadata) 

df_plot %>% 
  ggplot(aes(y = !!sym(key_gene), 
             x = ID)) + 
  geom_violin() + 
  # geom_boxplot(alpha = 0.5, outliers = FALSE) + 
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(
    title = glue('Per-condition expression level of {key_gene}/{names(key_gene)}')
  )

ggsave(glue('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/per_condition_expression_levels_{names(key_gene)}.pdf'),
       width = 10,
       height = 6)
  
# 🐸	Barplot of number of immunomodulatory genes significantly upregulated per condition.

DEG_up_immunomod_genes <- purrr::map_dfr(names(res), function(coef) {
  res[[coef]] %>%
    mutate(
      gene = sub("\\..*$", "", gene),  # strip version
      condition = coef
    ) %>%
    filter(gene %in% immunomod_genes)
    # filter(gene %in% immunomod_genes, padj < 0.05, log2FoldChange > 1) # Only significant up-regulated hits
})


# Count number of upregulated genes per condition
counts_per_condition <- DEG_up_immunomod_genes %>%
  group_by(condition) %>%
  summarise(n_up = n())

# Plot as barplot
ggplot(counts_per_condition, 
       aes(x = condition, 
           y = n_up)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Number of immunomodulatory genes significantly upregulated per condition",
    x = "Condition",
    y = "Number of upregulated genes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/immunomodulatory_upregulated_per_condition.pdf',
       width = 12,
       height = 6)
  

# Table of immunomodulatory DEGs per condition.
DEG_up_immunomod_genes_list <- DEG_up_immunomod_genes %>%
  group_split(condition, .keep = TRUE) %>%
  setNames(unique(DEG_up_immunomod_genes$condition)) %>%
  as.list()  # <-- converts it to a plain list

# Save output
# Path to output file
out_file <- "06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/out/DEG_up_immunomod.xlsx"

# If file exists, delete it first (optional)
if (file.exists(out_file)) file.remove(out_file)

# Loop over your list of contrasts
# Set Java memory BEFORE loading xlsx
# options(java.parameters = "-Xmx4g")
# Sys.setenv(JAVA_HOME = "/Library/Java/JavaVirtualMachines/temurin-17.jdk/Contents/Home")
library(rJava)
library(xlsx)
for (name in names(DEG_up_immunomod_genes_list)){
  # print((name))
  write.xlsx(DEG_up_immunomod_genes_list[[name]] %>% as.data.frame(), # Needs to be dataframe, cannot handle tibble
             file = out_file,
             sheetName = name,
             append = TRUE,
             row.names = FALSE)
}


  