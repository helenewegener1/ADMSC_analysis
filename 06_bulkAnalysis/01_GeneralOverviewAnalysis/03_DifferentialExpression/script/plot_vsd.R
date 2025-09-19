setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

library(DESeq2)
library(apeglm)
library(glue)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

dds <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds.rds')

rowMeans(counts(dds)) %>% hist(breaks = 50)

vsd_blind <- vst(dds, blind = TRUE)
vsd <- vst(dds, blind = FALSE)

vsd_mat_blind <- assay(vsd_blind)
vsd_mat <- assay(vsd)  # rows = genes, columns = samples

################################### PLOT VSD ################################### 

density_plot_vsd <- function(df){
  
  p <- df %>% 
    as.data.frame() %>% 
    rownames_to_column('genes') %>% 
    pivot_longer(cols = -genes) %>% 
    ggplot(aes(x = value, 
               color = name)) + 
    geom_density() +
    theme_bw() + 
    theme(legend.position = 'none')
  
  print(p)
  
}

boxplot_vsd <- function(df){
  
  p <- df %>% 
    as.data.frame() %>% 
    rownames_to_column('genes') %>% 
    pivot_longer(cols = -genes) %>% 
    ggplot(aes(y = value, 
               x = name)) + 
    geom_boxplot() +
    theme_bw() + 
    theme(legend.position = 'none')
  
  print(p)
  
}

################################################################################ 

density_plot_vsd(vsd_mat_blind)
density_plot_vsd(vsd_mat)


# boxplot_vsd(vsd_mat_blind)
# boxplot_vsd(vsd_mat)


################################## Filtering ###################################

# Filter dds 
smallestGroupSize <- 5 # each treatment has 5 donors
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds_filtered <- dds[keep,]

# Run vst
vsd_filtered_blind <- vst(dds_filtered, blind = TRUE)
vsd_filtered <- vst(dds_filtered, blind = FALSE)

vsd_filtered_mat_blind <- assay(vsd_filtered_blind)
vsd_filtered_mat <- assay(vsd_filtered)  # rows = genes, columns = samples

# Plot
density_plot_vsd(vsd_filtered_mat_blind)
density_plot_vsd(vsd_filtered_mat)

boxplot_vsd(vsd_filtered_mat_blind)
boxplot_vsd(vsd_filtered_mat)


