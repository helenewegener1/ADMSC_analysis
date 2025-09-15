
venn_diagram_synergy <- function(comparison, name_comparison){
  
  DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  
  if (length(DEGs_A) == 0 & length(DEGs_B) == 0){
    return('No DEGs to compare')
  }
  
  if (length(comparison) == 3){
    DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
    
    pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/venn_diagram_{name_comparison}.pdf'))
    grid.newpage() 
    venn.plot <- venn.diagram(
      x = list(A = DEGs_A, B = DEGs_B, C = DEGs_C),
      filename = NULL,  # NULL to plot in R
      alpha = 0.5,
      fill = c('darkred', 'darkblue', 'darkgreen'),
      category.names = names(comparison),
      cex = 1.5,
      cat.cex = 1,
      main = glue("Venn Diagram of DEGs between {name_comparison} VS ..."),
      disable.logging = TRUE
    )
    grid.draw(venn.plot)
    dev.off()
    
  } else if (length(comparison) == 2){
    
    pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/venn_diagram_{name_comparison}.pdf'))
    grid.newpage() 
    venn.plot <- venn.diagram(
      x = list(A = DEGs_A, B = DEGs_B),
      filename = NULL,  # NULL to plot in R
      alpha = 0.5,
      fill = c('darkred', 'darkblue'),
      category.names = names(comparison),
      cex = 1.5,
      cat.cex = 1,
      main = glue("Venn Diagram of DEGs between {name_comparison} VS ..."),
      disable.logging = TRUE
    )
    grid.draw(venn.plot)
    dev.off()
    
  } 
  
}


heatmap_synergy <- function(comparison){
  
  DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  
  synergy_genes <- union(DEGs_A, DEGs_B)
  
  if (length(comparison) == 3){
    DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
    synergy_genes <- union(synergy_genes, DEGs_C)
  } 
  
  # Subset vsd for synergy genes
  vsd_mat_synergy <- vsd_mat[synergy_genes, ]
  
  # Assess samples of conditions 
  IDs <- names(comparison) %>% str_remove_all('\\...')
  IDs <- c(IDs, paste(IDs, collapse = "_"))
  samples_of_condition <- meta_data %>% as.data.frame() %>% filter(ID %in% IDs) %>% arrange(ID)
  
  # Subset mat_top to samples of condition 
  vsd_mat_synergy_subset <- vsd_mat_synergy[, colnames(vsd_mat_synergy) %in% rownames(samples_of_condition)] 
  vsd_mat_synergy_subset <- vsd_mat_synergy_subset[, rownames(samples_of_condition)]
  
  pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/heatmap_syngergy_{IDs[length(IDs)]}.pdf'))
  pheatmap(
    vsd_mat_synergy_subset, 
    annotation_col = samples_of_condition %>% select(ID),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    main = glue("Synergy DEGs in {IDs[length(IDs)]}")
  )
  dev.off()
  
}

clustered_heatmap_synergy <- function(comparison){
  
  DEGs_A <- comparison[[1]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  DEGs_B <- comparison[[2]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
  
  synergy_genes <- union(DEGs_A, DEGs_B)
  
  if (length(comparison) == 3){
    DEGs_C <- comparison[[3]] %>% filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% select(gene) %>% unlist()
    synergy_genes <- union(synergy_genes, DEGs_C)
  } 
  
  # Subset vsd for synergy genes
  vsd_mat_synergy <- vsd_mat[synergy_genes, ]
  
  # Order samples 
  meta_data_order <- meta_data %>% as.data.frame() %>% arrange(ID)
  vsd_mat_synergy <- vsd_mat_synergy[, rownames(meta_data_order)]
  
  name <- names(comparison) %>% str_remove_all('\\...') %>% paste(collapse = "_")
  
  pdf(glue('06_bulkAnalysis/02_CytokineSynergyAnalysis/plot/clustered_heatmap_syngergy_{name}.pdf'))
  pheatmap(
    vsd_mat_synergy,  
    annotation_col = meta_data_order,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    main = glue("Synergy DEGs in {name} (all conditions)")
  )
  dev.off()
  
}
