#!/bin/bash

## Step 08 (Illumina): adapter and quality trimming with Trimmomatic, then FastQC and MultiQC
## Script originally written by Bernadette Johnson
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally.
#
# Usage (from the project folder; Trimmomatic and FastQC are called directly, no env):
#   bash 01_build_longread_transcriptome/08_illumina_trimmomatic.sh
# (MultiQC runs from its own env, multiqc_env)
#
# Reads:  data/01_build_longread_transcriptome/07_cutadapt/<sample>_R{1,2}.fastq.gz
#         data/reference/TruSeq3-PE-2.fa
# Writes: data/01_build_longread_transcriptome/08_trimmomatic/
#           <sample>_R{1,2}_paired.fastq.gz   trimmed pairs, input to steps 09 and 10
#           unpaired/                         reads whose mate was dropped
#           TrimmomaticOutput.txt             Trimmomatic log for all samples
#           fastqc/                           FastQC reports + multiqc_report.html
#
# A sample is skipped if its _R2_paired.fastq.gz already exists. Delete a
# sample's files to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/07_cutadapt
OUT=$DATA/08_trimmomatic
ADAPTERS=$ROOT/data/reference/TruSeq3-PE-2.fa
TRIMMOMATIC=~/Documents/programs/Trimmomatic-0.39/trimmomatic-0.39.jar
THREADS=16          # speed only, does not change results

[[ -f "$ADAPTERS" ]] || { echo "ERROR: adapter file $ADAPTERS not found" >&2; exit 1; }
[[ -f "$TRIMMOMATIC" ]] || { echo "ERROR: Trimmomatic not found at $TRIMMOMATIC" >&2; exit 1; }

shopt -s nullglob
r1_files=("$IN"/*_R1.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1.fastq.gz files in $IN. Run step 07 first." >&2; exit 1; }

mkdir -p "$OUT/unpaired" "$OUT/fastqc"
echo "Starting Trimmomatic..."
date

for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1.fastq.gz)

    if [[ -f "$OUT/${sample}_R2_paired.fastq.gz" ]]; then
        echo "$sample already done, skipping"
        continue
    fi

    echo "Processing: $sample"
    echo "=== $sample" >> "$OUT/TrimmomaticOutput.txt"
    # Crop the last 5 bp (CROP:145) because cutadapt poly-G trimming introduces a
    # bias against G nucleotides at this end
    java -jar "$TRIMMOMATIC" PE \
        -threads "$THREADS" \
        "$IN/${sample}_R1.fastq.gz" \
        "$IN/${sample}_R2.fastq.gz" \
        "$OUT/${sample}_R1_paired.fastq.gz" \
        "$OUT/unpaired/${sample}_R1_unpaired.fastq.gz" \
        "$OUT/${sample}_R2_paired.fastq.gz" \
        "$OUT/unpaired/${sample}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:2:keepBothReads \
        LEADING:3 \
        SLIDINGWINDOW:4:15 \
        CROP:145 \
        MINLEN:25 \
        >> "$OUT/TrimmomaticOutput.txt" 2>&1
done

echo "Completed Trimmomatic."
date

# FastQC
for fq in "$OUT"/*_paired.fastq.gz; do
    report=$OUT/fastqc/$(basename "$fq" .fastq.gz)_fastqc.html
    [[ -f "$report" ]] && continue
    echo "Running FastQC on $(basename "$fq")..."
    fastqc "$fq" -o "$OUT/fastqc"
done

# MultiQC
conda run -n multiqc_env multiqc -f "$OUT/fastqc" -o "$OUT/fastqc"

echo "All done! Thank you, please come again."
date
