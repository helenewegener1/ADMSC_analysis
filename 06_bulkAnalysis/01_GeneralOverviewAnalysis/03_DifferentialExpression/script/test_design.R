setwd('~/Documents/projects/project_ADMSC/analysis/')

# 3.	Differential expression
# o	For each cytokine condition vs Control (within each donor, then combined across donors if appropriate).
# HW: within each donor does not make sense since you are comparing 1 against 1 and then you have no dispersion. 
# o	For combined cytokine stimulations vs their respective single cytokines (to assess synergy).

library(DESeq2)
library(apeglm)
library(glue)
library(stringr)
library(ggplot2)

dds_1 <- DESeqDataSetFromTximport(txi, colData = metadata_clean, design = ~ ID )
dds_2 <- DESeqDataSetFromTximport(txi, colData = metadata_clean, design = ~ ID + Donor) # BEST BY FAR
dds_3 <- DESeqDataSetFromTximport(txi, colData = metadata_clean, design = ~ ID + Batch)
# dds_4 <- DESeqDataSetFromTximport(txi, colData = metadata_clean, design = ~ ID + Donor + Batch) # Model matrix not full rank

smallestGroupSize <- 5 # each treatment has 5 donors
keep <- rowSums(counts(dds_1) >= 10) >= smallestGroupSize
table(keep)

dds_1 <- dds_1[keep,]
dds_2 <- dds_2[keep,]
dds_3 <- dds_3[keep,]

dds_1 <- DESeq(dds_1)
res_1 <- results(dds_1)

dds_2 <- DESeq(dds_2)
res_2 <- results(dds_2)

dds_3 <- DESeq(dds_3)
res_3 <- results(dds_3)


# Contrasts
resultsNames(dds_1)

# ID_IFNy_vs_Ctrl
res_1_A <- results(dds_1, name = "ID_IFNy_vs_Ctrl")
res_1_A %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_2_A <- results(dds_2, name = "ID_IFNy_vs_Ctrl")
res_2_A %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_3_A <- results(dds_3, name = "ID_IFNy_vs_Ctrl")
res_3_A %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

# ID_IFNy_TNFa_IL17a_vs_Ctrl
res_1_B <- results(dds_1, name = "ID_IFNy_TNFa_IL17a_vs_Ctrl")
res_1_B %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_2_B <- results(dds_2, name = "ID_IFNy_TNFa_IL17a_vs_Ctrl")
res_2_B %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_3_B <- results(dds_3, name = "ID_IFNy_TNFa_IL17a_vs_Ctrl")
res_3_B %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

# IFNy_TNFa_vs_TNFa
res_1_C <- results(dds_1, contrast = c('ID', 'IFNy_TNFa', 'TNFa'))
res_1_C %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_2_C <- results(dds_2, contrast = c('ID', 'IFNy_TNFa', 'TNFa'))
res_2_C %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

res_3_C <- results(dds_3, contrast = c('ID', 'IFNy_TNFa', 'TNFa'))
res_3_C %>% ggplot(aes(pvalue)) + geom_histogram(binwidth = 0.05, boundary = 0, closed = "right")

