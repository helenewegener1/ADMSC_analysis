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

# Load data
dds <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds.rds')
# dds_synergy <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_synergy.rds')
res <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res.rds')
# res_synergy <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res_synergy.rds')
# resLFC <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/resLFC.rds')
meta_data <- dds@colData

source('./06_bulkAnalysis/02_CytokineSynergyAnalysis/script/functions.R')

# 🐝 Venn diagrams of DEGs between single and combined stimulation.

# comparison 1:	IFNγ + IL17A vs IFNγ and vs IL17A.
comparison1 <- list(res$IFNy_IL17a_vs_IFNy, res$IFNy_IL17a_vs_IL17a)
names(comparison1) <- c("...IFNy", "...IL17a")
name_comparison1 <- "IFNy_IL17a"

# comparison 2:	IFNγ + TNFα vs IFNγ and vs TNFα.
comparison2 <- list(res$IFNy_TNFa_vs_IFNy, res$IFNy_TNFa_vs_TNFa)
names(comparison2) <- c("...IFNy", "...TNFa")
name_comparison2 <- "IFNy_TNFa"

# comparison 3:	TNFα + IL17A vs TNFα and vs IL17A.
comparison3 <- list(res$TNFa_IL17a_vs_TNFa, res$TNFa_IL17a_vs_IL17a)
names(comparison3) <- c("...TNFa", "...IL17a")
name_comparison3 <- "TNFa_IL17a"

# comparison 4:	IFNγ + TNFα + IL17A vs each single cytokine and vs each pair.
comparison4 <- list(res$IFNy_TNFa_IL17a_vs_IFNy, res$IFNy_TNFa_IL17a_vs_IL17a, res$IFNy_TNFa_IL17a_vs_TNFa)
names(comparison4) <- c("...IFNy", "...IL17a", "...TNFa")
name_comparison4 <- "IFNy_TNFa_IL17a"

venn_diagram_synergy(comparison1, name_comparison1)
venn_diagram_synergy(comparison2, name_comparison2)
venn_diagram_synergy(comparison3, name_comparison3)
venn_diagram_synergy(comparison4, name_comparison4)

# Run with LFC-shrinken - shows no DEG's 
# venn_diagram_synergy(comparison = list("...IFNy" = resLFC$IFNy_IL17a_vs_IFNy, "...IL17a" = resLFC$IFNy_IL17a_vs_IL17a), 
#                      name_comparison = 'synergy_IFNy_IL17a')
# 
# venn_diagram_synergy(comparison = list("...IFNy" = resLFC$IFNy_TNFa_vs_IFNy, "...TNFa" = resLFC$IFNy_TNFa_vs_TNFa), 
#                      name_comparison = 'synergy_IFNy_TNFa')
# 
# venn_diagram_synergy(comparison = list("...TNFa" = resLFC$TNFa_IL17a_vs_TNFa, "...IL17a" = resLFC$TNFa_IL17a_vs_IL17a), 
#                      name_comparison = 'synergy_TNFa_IL17a')



# 🐝 Heatmaps of synergy-specific genes.

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Comparisons of interest:
# comparison 1:	IFNγ + IL17A vs IFNγ and vs IL17A.
# comparison 2:	IFNγ + TNFα vs IFNγ and vs TNFα.
# comparison 3:	TNFα + IL17A vs TNFα and vs IL17A.
# comparison 4:	IFNγ + TNFα + IL17A vs each single cytokine and vs each pair.

heatmap_synergy(comparison = comparison1)
heatmap_synergy(comparison = comparison2)
heatmap_synergy(comparison = comparison3)
heatmap_synergy(comparison = comparison4)


# 🐝 Clustered heatmap showing expression patterns across all conditions.

clustered_heatmap_synergy(comparison = comparison1)
clustered_heatmap_synergy(comparison = comparison2)
clustered_heatmap_synergy(comparison = comparison3)
clustered_heatmap_synergy(comparison = comparison4)
