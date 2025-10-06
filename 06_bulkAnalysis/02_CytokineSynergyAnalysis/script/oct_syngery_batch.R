setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

# Comparisons of interest:
# comparison 1:	IFNγ + IL17A vs IFNγ and vs IL17A.
# comparison 2:	IFNγ + TNFα vs IFNγ and vs TNFα.
# comparison 3:	TNFα + IL17A vs TNFα and vs IL17A.
# comparison 4:	IFNγ + TNFα + IL17A vs each single cytokine and vs each pair.

# TODO: Table of synergy-specific genes. (working on it in functions )

library(dplyr)
library(pheatmap)
library(glue)
library(DESeq2)
library(stringr)
library(ggplot2)
library(VennDiagram)
library(biomaRt)

# batch_nr <- 1
# batch_nr <- 2

dds <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))
meta_data <- dds@colData

source('./06_bulkAnalysis/02_CytokineSynergyAnalysis/script/functions.R')

# Make sure Control is the reference level
dds$ID <- relevel(dds$ID, ref = "Ctrl")

# Re-run DESeq2
dds <- DESeq(dds)

# Check names
resultsNames(dds)

oct_synergy_heatmap("IFNy", "IL17a")
oct_synergy_heatmap("IFNy", "TNFa")
oct_synergy_heatmap("IFNy", "TNFa")
oct_synergy_heatmap("TNFa", "IL17a")
oct_synergy_heatmap("IFNy", "TNFa", "IL17a")

# # Define synergy to analyze 
# cytokine1 <- "IFNy"
# cytokine2 <- "IL17a"
# cytokine3 <- NULL
# cytokine3 <- "TNFa"
# 
# if (!is.null(cytokine3)){
#   cytokine_combo <- paste(cytokine1, cytokine2, cytokine3, sep = "_")
# } else {
#   cytokine_combo <- paste(cytokine1, cytokine2, sep = "_")
# }
# 
# # Build contrast vector
# all_coefs <- resultsNames(dds)
# contrast_vec <- numeric(length(all_coefs))
# names(contrast_vec) <- all_coefs
# 
# contrast_vec[glue("ID_{cytokine_combo}_vs_Ctrl")] <- 1
# contrast_vec[glue("ID_{cytokine1}_vs_Ctrl")] <- -1
# contrast_vec[glue("ID_{cytokine2}_vs_Ctrl")] <- -1
# 
# if (!is.null(cytokine3)){
#   contrast_vec[glue("ID_{cytokine3}_vs_Ctrl")] <- -1
# } 
# 
# # Test synergy
# synergy <- results(dds, contrast = contrast_vec)
# 
# # Significant synergistic genes
# sig_synergy <- subset(synergy, padj < 0.05)
# positive_synergy <- subset(sig_synergy, log2FoldChange > 1)  # Amplification
# negative_synergy <- subset(sig_synergy, log2FoldChange < 1)  # Dampening
# 
# summary(synergy)
# 
# # Extract list of synergy genes
# synergy_genes <- c(rownames(positive_synergy), rownames(negative_synergy))
# 
# # Heatmap 
# vsd <- vst(dds, blind = FALSE)
# vsd_mat <- assay(vsd)
# 
# # Subset vsd for synergy genes
# vsd_mat_synergy <- vsd_mat[synergy_genes, ]
# 
# # Order samples
# levels <- str_split(cytokine_combo, "_") %>% unlist() %>% c(cytokine_combo, "Ctrl")
# 
# meta_data_order <- meta_data %>% as.data.frame() %>% 
#   filter(ID %in% levels) %>% select(ID) %>% 
#   mutate(ID = factor(ID, levels = levels)) %>% 
#   arrange(ID)
# 
# vsd_mat_synergy <- vsd_mat_synergy[, rownames(meta_data_order)]
# 
# pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/oct_synergy/heatmap_{cytokine_combo}_synergy.pdf'), height = 12)
# pheatmap(
#   vsd_mat_synergy,  
#   annotation_col = meta_data_order,
#   cluster_cols = FALSE,
#   cluster_rows = FALSE,
#   show_colnames = FALSE
#   # main = glue("Synergy DEGs in {name} (all conditions)")
# )
# dev.off()



