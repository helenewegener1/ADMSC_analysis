# 
# venn_diagram_synergy <- function(comparison, name_comparison){
#   
#   DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
#   DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
#   
#   if (length(DEGs_A) == 0 & length(DEGs_B) == 0){
#     return('No DEGs to compare')
#   }
#   
#   if (length(comparison) == 3){
#     DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
#     
#     pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/venn_diagram_{name_comparison}.pdf'))
#     grid.newpage() 
#     venn.plot <- venn.diagram(
#       x = list(A = DEGs_A, B = DEGs_B, C = DEGs_C),
#       filename = NULL,  # NULL to plot in R
#       alpha = 0.5,
#       fill = c('darkred', 'darkblue', 'darkgreen'),
#       category.names = names(comparison),
#       cex = 1.5,
#       cat.cex = 1,
#       main = glue("Venn Diagram of DEGs between {name_comparison} VS ..."),
#       disable.logging = TRUE
#     )
#     grid.draw(venn.plot)
#     dev.off()
#     
#   } else if (length(comparison) == 2){
#     
#     pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/venn_diagram_{name_comparison}.pdf'))
#     grid.newpage() 
#     venn.plot <- venn.diagram(
#       x = list(A = DEGs_A, B = DEGs_B),
#       filename = NULL,  # NULL to plot in R
#       alpha = 0.5,
#       fill = c('darkred', 'darkblue'),
#       category.names = names(comparison),
#       cex = 1.5,
#       cat.cex = 1,
#       main = glue("Venn Diagram of DEGs between {name_comparison} VS ..."),
#       disable.logging = TRUE
#     )
#     grid.draw(venn.plot)
#     dev.off()
#     
#   } 
#   
# }
# 
# 
# heatmap_synergy <- function(comparison){
#   
#   DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   
#   synergy_genes <- union(DEGs_A, DEGs_B)
#   
#   if (length(comparison) == 3){
#     DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#     synergy_genes <- union(synergy_genes, DEGs_C)
#   } 
#   
#   # Subset vsd for synergy genes
#   vsd_mat_synergy <- vsd_mat[synergy_genes, ]
#   
#   # Assess samples of conditions 
#   IDs <- names(comparison) %>% str_remove_all('\\...')
#   IDs <- c(IDs, paste(IDs, collapse = "_"))
#   samples_of_condition <- meta_data %>% as.data.frame() %>% filter(ID %in% IDs) %>% arrange(ID)
#   
#   # Subset mat_top to samples of condition 
#   vsd_mat_synergy_subset <- vsd_mat_synergy[, colnames(vsd_mat_synergy) %in% rownames(samples_of_condition)] 
#   vsd_mat_synergy_subset <- vsd_mat_synergy_subset[, rownames(samples_of_condition)]
#   
#   pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/heatmap_syngergy_{IDs[length(IDs)]}.pdf'))
#   pheatmap(
#     vsd_mat_synergy_subset, 
#     annotation_col = samples_of_condition %>% select(ID),
#     cluster_cols = FALSE,
#     cluster_rows = FALSE,
#     main = glue("Synergy DEGs in {IDs[length(IDs)]}")
#   )
#   dev.off()
#   
# }
# 
# clustered_heatmap_synergy <- function(comparison){
#   
#   DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   
#   synergy_genes <- union(DEGs_A, DEGs_B)
#   
#   if (length(comparison) == 3){
#     DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#     synergy_genes <- union(synergy_genes, DEGs_C)
#   } 
#   
#   # Subset vsd for synergy genes
#   vsd_mat_synergy <- vsd_mat[synergy_genes, ]
#   
#   # Order samples 
#   meta_data_order <- meta_data %>% as.data.frame() %>% arrange(ID)
#   vsd_mat_synergy <- vsd_mat_synergy[, rownames(meta_data_order)]
#   
#   name <- names(comparison) %>% str_remove_all('\\...') %>% paste(collapse = "_")
#   
#   pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/clustered_heatmap_syngergy_{name}.pdf'))
#   pheatmap(
#     vsd_mat_synergy,  
#     annotation_col = meta_data_order,
#     cluster_cols = TRUE,
#     cluster_rows = TRUE,
#     show_colnames = FALSE,
#     main = glue("Synergy DEGs in {name} (all conditions)")
#   )
#   dev.off()
#   
# }
# 
# 
# # TODO: FIGURE OUT HOW THIS SHOULD BE DONE
# table_synergy <- function(comparison, name_comparison){
#   
#   DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% dplyr::select(gene) %>% unlist()
#   
#   synergy_genes <- union(DEGs_A, DEGs_B)
#   
#   if (length(comparison) == 3){
#     DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
#     synergy_genes <- union(synergy_genes, DEGs_C)
#   } 
#   
#   # Assess samples of conditions 
#   IDs <- names(comparison) %>% str_remove_all('\\...')
#   IDs <- c(paste(IDs, collapse = "_"), IDs)
#   
#   compare1 <- glue('{IDs[1]}_vs_{IDs[2]}')
#   compare2 <- glue('{IDs[1]}_vs_{IDs[3]}')
#   
#   if (length(comparison) == 3){
#     compare3 <- glue('{IDs[1]}_vs_{IDs[3]}')
#   } 
#   
#   # Subset synergy genes in res
#   res[[compare1]] %>% filter(gene %in% synergy_genes) 
#   
# }


oct_synergy_heatmap <- function(cytokine1, cytokine2, cytokine3 = NULL) {
  
  if (!is.null(cytokine3)){
    cytokine_combo <- paste(cytokine1, cytokine2, cytokine3, sep = "_")
  } else {
    cytokine_combo <- paste(cytokine1, cytokine2, sep = "_")
  }
  
  # Build contrast vector
  all_coefs <- resultsNames(dds)
  contrast_vec <- numeric(length(all_coefs))
  names(contrast_vec) <- all_coefs
  
  contrast_vec[glue("ID_{cytokine_combo}_vs_Ctrl")] <- 1
  contrast_vec[glue("ID_{cytokine1}_vs_Ctrl")] <- -1
  contrast_vec[glue("ID_{cytokine2}_vs_Ctrl")] <- -1
  
  if (!is.null(cytokine3)){
    contrast_vec[glue("ID_{cytokine3}_vs_Ctrl")] <- -1
  } 
  
  # Test synergy
  synergy <- results(dds, contrast = contrast_vec)
  
  # Significant synergistic genes
  sig_synergy <- subset(synergy, padj < 0.05)
  positive_synergy <- subset(sig_synergy, log2FoldChange > 1) # Amplification
  negative_synergy <- subset(sig_synergy, log2FoldChange < 1)  # Dampening
  
  summary(synergy)
  
  positive_synergy_top50 <- positive_synergy %>% as.data.frame() %>% arrange(padj) %>% head(50)
  negative_synergy_top50 <- negative_synergy %>% as.data.frame() %>% arrange(padj) %>% head(50)
  
  # Extract list of synergy genes
  synergy_genes <- c(rownames(negative_synergy_top50), rownames(positive_synergy_top50))
  
  # Heatmap 
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  
  # Subset vsd for synergy genes
  vsd_mat_synergy <- vsd_mat[synergy_genes, ]
  
  # Order samples
  # levels <- str_split(cytokine_combo, "_") %>% unlist() %>% c(cytokine_combo, "Ctrl")
  
  meta_data_order <- meta_data %>% as.data.frame() %>% 
    arrange(ID) %>% 
    dplyr::select(ID)
  #   filter(ID %in% levels) %>% select(ID, Donor) %>% 
  #   mutate(ID = factor(ID, levels = levels)) %>%
  #   arrange(ID)
  
  
  # Connect to Ensembl
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get mapping
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = synergy_genes %>% str_split_i('\\.', 1),
    mart = mart
  )
  
  gene_map <- gene_map %>% mutate(hgnc_symbol = case_when(hgnc_symbol == "" ~ ensembl_gene_id, .default = hgnc_symbol))
  
  rownames(vsd_mat_synergy) <- rownames(vsd_mat_synergy) %>% str_split_i('\\.', 1)
  
  vsd_mat_synergy <- vsd_mat_synergy[gene_map$ensembl_gene_id, ]
  
  # check order 
  table(rownames(vsd_mat_synergy) == gene_map$ensembl_gene_id)
  
  # Rename genes 
  rownames(vsd_mat_synergy) <- gene_map$hgnc_symbol

  vsd_mat_synergy <- vsd_mat_synergy[, rownames(meta_data_order)]
  
  pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/oct_synergy/batch_{batch_nr}_heatmap_{cytokine_combo}_synergy.pdf'), height = 12)
  pheatmap(
    vsd_mat_synergy,  
    annotation_col = meta_data_order,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE
    # main = glue("Synergy DEGs in {name} (all conditions)")
  )
  dev.off()
  
}
