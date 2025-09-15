setwd('~/Documents/projects/project_ADMSC/ADMSC_analysis/')

library(rtracklayer) 

# Replace with the path to your downloaded GTF
gtf_file <- "../data/human_reference/gencode.v49.basic.annotation.gtf.gz"

# Import GTF
gtf <- import(gtf_file)

# Keep only transcript features
transcripts <- gtf[gtf$type == "transcript"]

# Create tx2gene table
tx2gene <- data.frame(
  transcript_id = transcripts$transcript_id,
  gene_id = transcripts$gene_id
)

# Save to tab-separated file
write.table(tx2gene, file="05_tx2gene/out/tx2gene.tsv", sep="\t", quote=FALSE, row.names=FALSE)

