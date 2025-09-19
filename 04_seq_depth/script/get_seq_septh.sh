#!/bin/bash

for f in /Users/srz223/Documents/projects/project_ADMSC/ADMSC_analysis/03_kallisto_quant/out/*/abundance.tsv; do
    sample=$(basename $(dirname "$f"))
    total=$(awk 'NR>1 {sum += $4} END {print sum}' "$f")
    echo -e "$sample\t$total"
done > ../out/sequencing_depths.txt

