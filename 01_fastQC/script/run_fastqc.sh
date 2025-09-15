#!/bin/bash

cd /Users/srz223/Documents/projects/project_ADMSC/data/RAW_DATA

samples=$(ls)

outdir="/Users/srz223/Documents/projects/project_ADMSC/ADMSC_analysis/01_fastQC/out"
mkdir -p "$outdir"

for sample in $(ls); do

    fastqc --threads 6 "${sample}/${sample}_1.fq.gz" "${sample}/${sample}_2.fq.gz" --outdir "$outdir"

done


# 
