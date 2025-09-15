#!/bin/bash

# Paths
RAW_DIR="/Users/srz223/Documents/projects/project_ADMSC/data/RAW_DATA"
OUT_DIR="../out"
IDX="/Users/srz223/Documents/projects/project_ADMSC/ADMSC_analysis/02_kallisto_index/out/kallisto_human49_index.idx"

mkdir -p "$OUT_DIR"

THREADS=6      # threads per kallisto job
MAX_JOBS=2     # number of parallel samples

# Function: wait until fewer than MAX_JOBS are running
wait_for_jobs() {
    while (( $(jobs -r | wc -l) >= MAX_JOBS )); do
        sleep 1
    done
}

for sample_dir in "$RAW_DIR"/*/; do
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"

    if [[ "$sample" == "4" ]]; then
        echo "Running sample 4 as single-end"
        kallisto quant -i "$IDX" \
                       -o "$OUT_DIR/$sample" \
                       -b 100 \
                       -t $THREADS \
                       --single -l 200 -s 20 \
                       "$sample_dir/${sample}_1.fq.gz" &
    elif [[ -f "$sample_dir/${sample}_1.fq.gz" && -f "$sample_dir/${sample}_2.fq.gz" ]]; then
        kallisto quant -i "$IDX" \
                       -o "$OUT_DIR/$sample" \
                       -b 100 \
                       -t $THREADS \
                       "$sample_dir/${sample}_1.fq.gz" "$sample_dir/${sample}_2.fq.gz" &
    else
        echo "No paired-end files found for $sample, skipping..."
    fi

    wait_for_jobs   # waits until fewer than MAX_JOBS are running
done

wait
echo "All samples processed"
