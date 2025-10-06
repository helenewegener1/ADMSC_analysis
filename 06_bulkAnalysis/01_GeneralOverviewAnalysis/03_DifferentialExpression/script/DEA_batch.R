setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

# 3.	Differential expression
# o	For each cytokine condition vs Control (within each donor, then combined across donors if appropriate).
# HW: within each donor does not make sense since you are comparing 1 against 1 and then you have no dispersion. 
# o	For combined cytokine stimulations vs their respective single cytokines (to assess synergy).

library(DESeq2)
library(apeglm)
library(glue)
library(stringr)
library(dplyr)
# detach("package:biomaRt", unload = TRUE, character.only = TRUE)

# batch_nr <- 1
# batch_nr <- 2

# Load DESeqDataSet object
dds <- readRDS(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_batch_{batch_nr}.rds'))

# Pre-filtering (todo or not?)
# smallestGroupSize <- table(metadata_clean$ID, metadata_clean$Donor) %>% min()
smallestGroupSize <- 5 # each treatment has 5 donors
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]

# Specify reference
dds$ID <- relevel(dds$ID, ref = "Ctrl")

################################################################################
####################### Differential expression analysis ####################### 
################################################################################

# Run DEA
dds <- DESeq(dds)
res <- results(dds)
res

# Contrasts
resultsNames(dds)

# 🐒 For each cytokine condition vs Control (within each donor, then combined across donors if appropriate).
# HW: within each donor does not make sense since you are comparing 1 against 1 and then you have no dispersion. 

# Prep output 
res <- list()
resLFC <- list()

# Assess all cytokine conditions
cytokine_conditions <- levels (dds$ID)[2:length(levels(dds$ID))] 

for (cytokine_condition in cytokine_conditions){
  
  name <- glue('{cytokine_condition}_vs_Ctrl')
  res_tmp <- results(dds, name = glue("ID_{cytokine_condition}_vs_Ctrl"))
  
  #Remove NA rows and filter by adjusted p value
  res_tmp <- res_tmp[!is.na(res_tmp$padj) ,]
  res_tmp <- res_tmp[order(res_tmp$padj), ]
  
  res[[name]] <- as.data.frame(res_tmp) %>% rownames_to_column('gene')
  
  # Log fold change shrinkage
  resLFC_tpm <- lfcShrink(dds, coef = glue("ID_{cytokine_condition}_vs_Ctrl"), type="apeglm")
  resLFC_tpm <- resLFC_tpm[!is.na(resLFC_tpm$padj) ,]
  resLFC_tpm <- resLFC_tpm[order(resLFC_tpm$padj), ]
  
  resLFC[[name]] <- as.data.frame(resLFC_tpm) %>% rownames_to_column('gene')
  
}

# 🐒	For combined cytokine stimulations vs their respective single cytokines (to assess synergy).

# Assess individual cytokines
cytokines <- c('IFNy', 'IL17a', 'TNFa')

for (cytokine in cytokines){
  for (cytokine_condition in cytokine_conditions){

    if ((str_detect(cytokine_condition, cytokine)) & (cytokine != cytokine_condition)){

      # Relevel
      dds$ID <- relevel(dds$ID, ref = cytokine)
      dds <- DESeq(dds)
      resultsNames(dds)

      # contrast = c(condition, level_of_interest, reference)
      name <- glue('{cytokine_condition}_vs_{cytokine}')

      res_tmp <- results(dds, contrast = c("ID", cytokine_condition, cytokine))

      #Remove NA rows and filter by adjusted p value
      res_tmp <- res_tmp[!is.na(res_tmp$padj) ,]
      res_tmp <- res_tmp[order(res_tmp$padj), ]

      res[[name]] <- as.data.frame(res_tmp) %>% rownames_to_column('gene')

      # Log fold change shrinkage
      resLFC_tpm <- lfcShrink(dds, coef = glue("ID_{cytokine_condition}_vs_{cytokine}"), type="apeglm")
      resLFC_tpm <- resLFC_tpm[!is.na(resLFC_tpm$padj) ,]
      resLFC_tpm <- resLFC_tpm[order(resLFC_tpm$padj), ]

      resLFC[[name]] <- as.data.frame(resLFC_tpm) %>% rownames_to_column('gene')



    }

  }
}

# Export res
saveRDS(res, glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res_batch_{batch_nr}.rds'))
saveRDS(resLFC, glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/resLFC_batch_{batch_nr}.rds'))


################################################################################
############################## P value histograms ############################## 
################################################################################

# Same plots somes with resLFC

comparisons <- names(res)

for(comparison in comparisons){
  
  res[[comparison]] %>% 
    ggplot(aes(pvalue)) + 
    geom_histogram(binwidth = 0.05, boundary = 0, closed = "right") + 
    theme_bw() + 
    labs(
      title = glue('P-value histogram of {comparison}')
    )
  
  ggsave(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/plot/batch_{batch_nr}_pvalue_histogram_{comparison}.pdf'), 
         width = 10,
         height = 6)
  
}











