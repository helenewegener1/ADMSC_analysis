#!/bin/bash

# Unfiltered FASTA file
kallisto index -i ../out/kallisto_human49_index.idx /Users/srz223/Documents/projects/project_ADMSC/data/human_reference/gencode.v49.transcripts.fa.gz

# Filtered FASTA file
kallisto index -i ../out/kallisto_human49_index_filtered.idx /Users/srz223/Documents/projects/project_ADMSC/data/human_reference/gencode.v49.transcripts.filtered.fa.gz


