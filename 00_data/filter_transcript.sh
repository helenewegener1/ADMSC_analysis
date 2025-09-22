cd /Users/srz223/Documents/projects/project_ADMSC/data

# Human transcriptome of all transcript sequences. 
# Release 49 (GRCh38.p14)
# website: https://www.gencodegenes.org/human/ 

# format: FASTA
# Transcript sequences - ALL
# Used to build index in 02_kallisto_index
# The FASTA provides the transcript sequences Kallisto uses for pseudoalignment.

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

# GFF file 
# Bacis annotation - CHR - GFT
# Used to build tx2gene in R
# The GTF defines which transcripts and genes exist in your analysis.
# The GTF contains a description of the coordinates of exons that make up each transcript 
# but it does not contain the transcript sequences themselves.

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz

# Here’s why filtering the FASTA is important for your Kallisto quantification:
#
# Consistency between annotation (GTF) and reference (FASTA):
# Your GTF defines which transcripts and genes exist in your analysis.
# Your FASTA provides the transcript sequences Kallisto uses for pseudoalignment.
# If the FASTA has transcripts not in the GTF (or vice versa), 
# Kallisto will quantify extra transcripts that have no annotation — leading to “unannotated” results or empty gene names downstream.
# 
# Avoids noisy counts:
# Including unannotated/patch/haplotype transcripts can spread reads across multiple nearly identical transcripts.
# That makes your counts smaller per transcript and can reduce your power to detect DEGs.
# 
# Improves interpretability:
# By filtering to basic transcripts, you’re focusing on well-supported, canonical transcripts.
# This avoids spurious low-quality transcripts and reduces false positives/negatives in DEG analysis.
# 
# Better comparability across samples & studies:
# Using GTF + matching FASTA from the same Gencode release ensures reproducibility and that your results will match published data.
# So yes — using a filtered FASTA matched to your chosen GTF is considered best practice for bulk RNA-seq (and single-cell too).

# Extract transcript IDs from the basic GTF
zgrep -v '^#' gencode.v49.basic.annotation.gtf.gz | \
    awk '$3=="transcript"{print $12}' | \
    sed 's/"//g; s/;//g' > basic_transcript_ids.txt

# Filter the FASTA to keep only those transcript IDs
gzcat gencode.v49.transcripts.fa.gz | \
awk 'BEGIN{
    while ((getline id < "basic_transcript_ids.txt") > 0) keep[id]=1
}
# When we see a header line
/^>/ {
    split($1, parts, "|")          # split at |
    sub(/^>/, "", parts[1])        # remove ">" from transcript ID
    keep_seq = (parts[1] in keep)  # check if transcript ID is in keep list
}
keep_seq' > gencode.v49.transcripts.filtered.fa

gzip gencode.v49.transcripts.filtered.fa 
