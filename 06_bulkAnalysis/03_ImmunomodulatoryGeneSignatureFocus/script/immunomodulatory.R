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

# batch_nr <- 1
batch_nr <- 2

# Load DESeqDataSet object
dds <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))
res <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res_batch_{batch_nr}.rds'))
metadata <- dds@colData %>% as.data.frame() %>% arrange(ID) 

# Anders + more added with AI
immunomod_genes_symbol <- c(
  # Your Initial List (Immunosuppression, MHC, Chemokines)
  "IDO1", "NOS2", "HHLA2", 
  "HLA-A", "HLA-B", "HLA-C", # HLA_A, _B, _C - class 1 - IFNy burde have effekt 
  "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", # HLADRPQ - class 2 - IFNy burde have effekt 
  "CXCL5", "CCL20", "CD274",
  
  # Core Cytokines & Chemokines (Pro/Anti-inflammatory)
  "IL2", "IL6", "IL10", "IFNG", "TNF", "TGFB1", "CCL2", "CXCL8", 
  
  # Core T-cell Checkpoints & Co-stimulation
  "PDCD1", "CTLA4", "CD28", "ICOS", 
  
  # Master Transcription Factors
  "NFKB1", "FOXP3", "TBX21", "GATA3", "RORC", 
  
  # Innate Sensors & Signaling
  "TLR4", "NLRP3", "MYD88"
)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(dds) %>% str_split_i('\\.', 1),
  mart = mart
)

# grep("CD274", gene_map$hgnc_symbol, value = TRUE, ignore.case = TRUE)

immunomod_genes <- gene_map %>% filter(hgnc_symbol %in% immunomod_genes_symbol) %>% dplyr::select(ensembl_gene_id) %>% unlist()
names(immunomod_genes) <- immunomod_genes_symbol

# 🐸	Heatmap of immunomodulatory gene expression across all samples.

# Get normalized expression
vsd <- vst(dds, blind = FALSE)
rownames(vsd) <- sub("\\..*$", "", rownames(vsd)) # Remove version numbers (everything after the dot)
mat <- assay(vsd)[rownames(vsd) %in% immunomod_genes, rownames(metadata)]

# check order 
table(rownames(mat) == immunomod_genes)

# Rename genes 
rownames(mat) <- names(immunomod_genes)

# Optionally order by gene family or variance
# mat <- mat[order(rowMeans(mat), decreasing = TRUE), ]

pdf(glue('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/batch_{batch_nr}_heatmap_immunomodulatory_all_samples.pdf'))

pheatmap(mat,
         annotation_col = metadata %>% dplyr::select(ID, Donor),
         scale = "row",
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_colnames = TRUE, 
         angle_col = 45,
         fontsize_row = 6,
         main = "Immunomodulatory Gene Expression")

dev.off()

# 🐸	Dot plot or boxplot of selected key genes showing per-condition expression levels.
key_genes <- list("IDO1", "NOS2", "HLA-A", "HLA-B", "HLA-C", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3")

df_plot <- mat %>% t() %>% cbind(metadata) 

for (key_gene in key_genes){
  
  # key_gene <- key_genes[[1]]
    
  p <- df_plot %>% 
    ggplot(aes(y = !!sym(key_gene), 
               x = ID)) + 
    stat_summary(fun = "mean", 
                 geom = "crossbar", 
                 width = 0.4,       
                 linewidth = 0.2) +
    geom_jitter(aes(color = Donor), width = 0.1, alpha = 0.5, size = 2) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(
      title = glue('Per-condition expression level of {key_gene}')
    )
  
  ggsave(glue('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/batch_{batch_nr}_per_condition_expression_levels_{key_gene}.pdf'),
         plot = p,
         width = 10,
         height = 6)
  
}


  
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

# ggsave(glue('06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/plot/batch_{batch_nr}_immunomodulatory_upregulated_per_condition.pdf'),
#        width = 12,
#        height = 6)
  

# Table of immunomodulatory DEGs per condition.
DEG_up_immunomod_genes_list <- DEG_up_immunomod_genes %>%
  group_split(condition, .keep = TRUE) %>%
  setNames(unique(DEG_up_immunomod_genes$condition)) %>%
  as.list()  # <-- converts it to a plain list

# Save output
# Path to output file
out_file <- glue("06_bulkAnalysis/03_ImmunomodulatoryGeneSignatureFocus/out/batch_{batch_nr}_DEG_up_immunomod.xlsx")

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


  