setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

# 3.	Differential expression
# o	For each cytokine condition vs Control (within each donor, then combined across donors if appropriate).
# HW: within each donor does not make sense since you are comparing 1 against 1 and then you have no dispersion. 
# o	For combined cytokine stimulations vs their respective single cytokines (to assess synergy).

library(DESeq2)
library(apeglm)
library(glue)
library(stringr)

# Load DESeqDataSet object
dds_synergy <- readRDS('06_bulkAnalysis/01_GeneralOverviewAnalysis/02_GlobalExpressionPatterns/out/dds_synergy.rds')

# Pre-filtering (todo or not?)
smallestGroupSize <- 5 # each treatment has 5 donors
keep <- rowSums(counts(dds_synergy) >= 10) >= smallestGroupSize
table(keep)
dds_synergy <- dds_synergy[keep,]

################################################################################
####################### Differential expression analysis ####################### 
################################################################################

# Run DEA
dds_synergy <- DESeq(dds_synergy)
res <- results(dds_synergy)
res

# Contrasts
resultsNames(dds_synergy)

# Main effects: Additive effect of each single cytokine:
# IFNy → Effect of IFNy alone, relative to baseline (Ctrl), ignoring the other cytokines.
# TNFa → Effect of TNFa alone, relative to baseline.
# IL17a → Effect of IL17a alone, relative to baseline.

# Two-way interactions: 
# Positive interaction coefficient → synergistic effect (greater than sum of singles).
# Negative → antagonistic (less than sum of singles).
# IFNy.TNFa → Additional effect of having both IFNy and TNFa together, beyond what would be expected if you just added the single effects.
# IFNy.IL17a → Additional effect of IFNy + IL17a together.
# TNFa.IL17a → Additional effect of TNFa + IL17a together.

# Three-way interaction
# This tells you if the triple combination has a unique effect that can’t be explained by any of the singles or pairs.
# IFNy.TNFa.IL17a → Extra effect of having all three cytokines together, beyond the sum of the single and pairwise effects.


# 🐒	For combined cytokine stimulations vs their respective single cytokines (to assess synergy).

# Prep output 
res <- list()
resLFC <- list()

# Assess all cytokine conditions
coefficients <- c("IFNy", "TNFa", "IL17a", 
                  "IFNy.TNFa", "IFNy.IL17a", "TNFa.IL17a",
                  "IFNy.TNFa.IL17a")

for (coefficient in coefficients){
  
  res_tmp <- results(dds_synergy, name = coefficient)
  
  #Remove NA rows and filter by adjusted p value
  res_tmp <- res_tmp[order(res_tmp$padj), ]
  res_tmp <- res_tmp[!is.na(res_tmp$padj) ,]
  
  res[[coefficient]] <- as.data.frame(res_tmp) %>% rownames_to_column('gene')
  
  # Log fold change shrinkage
  resLFC_tpm <- lfcShrink(dds_synergy, coef = coefficient, type="apeglm")
  resLFC_tpm <- resLFC_tpm[!is.na(resLFC_tpm$padj) ,]
  resLFC_tpm <- resLFC_tpm[order(resLFC_tpm$padj), ]
  
  resLFC[[coefficient]] <- as.data.frame(resLFC_tpm) %>% rownames_to_column('gene')
  
}

# Export res
saveRDS(res, '06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/res_synergy.rds')
saveRDS(resLFC, '06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/out/resLFC_synergy.rds')

################################################################################
############################## P value histograms ############################## 
################################################################################

# Same plots somes with resLFC

coefficients <- c("IFNy", "TNFa", "IL17a", 
                  "IFNy.TNFa", "IFNy.IL17a", "TNFa.IL17a",
                  "IFNy.TNFa.IL17a")

for (coefficient in coefficients){
  
  res[[coefficient]] %>% 
    ggplot(aes(pvalue)) + 
    geom_histogram(binwidth = 0.05, boundary = 0, closed = "right") + 
    theme_bw() + 
    labs(
      title = glue('P-value histogram of {coefficient}')
    )
  
  ggsave(glue('06_bulkAnalysis/01_GeneralOverviewAnalysis/03_DifferentialExpression/plot/synergy_pvalue_histogram_{coefficient}.pdf'), 
         width = 10,
         height = 6)
  
}








