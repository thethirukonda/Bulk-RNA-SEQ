#!/bin/bash
# Go to your aligned reads folder
cd /Users/aakash/Downloads/rnaseq/bam_files

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time
    echo "Processing $bam ..."
    featureCounts -s 2 -a /Users/aakash/Downloads/rnaseq/Homo_sapiens.GRCh38.115.gtf \
        -o /Users/aakash/Downloads/rnaseq/quants/${bam%.bam}_featurecounts.txt \
        "$bam"
    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes
    echo "Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
