#!/bin/bash

## Step 07 (Illumina): poly-G trimming with cutadapt, then FastQC and MultiQC
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally.
#
# Usage (from the project folder):
#   conda activate cutadapt
#   bash 01_build_longread_transcriptome/07_illumina_cutadapt.sh
# (MultiQC runs from its own env, multiqc_env)
#
# Reads:  data/01_build_longread_transcriptome/00_rawillumina/<sample>_R{1,2}.fastq.gz
# Writes: data/01_build_longread_transcriptome/07_cutadapt/
#           <sample>_R{1,2}.fastq.gz   trimmed reads, input to step 08
#           fastqc/                    FastQC reports + multiqc_report.html
#
# A sample is skipped if its trimmed _R2.fastq.gz already exists. Delete a
# sample's files to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/00_rawillumina
OUT=$DATA/07_cutadapt
THREADS=16          # speed only, does not change results

shopt -s nullglob
r1_files=("$IN"/*_R1.fastq.gz)
r2_files=("$IN"/*_R2.fastq.gz)
echo "Number of forward reads: ${#r1_files[@]}"
echo "Number of reverse reads: ${#r2_files[@]}"
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1.fastq.gz files in $IN" >&2; exit 1; }
[[ ${#r1_files[@]} -eq ${#r2_files[@]} ]] || { echo "ERROR: the number of forward reads does not match the number of reverse reads" >&2; exit 1; }

mkdir -p "$OUT/fastqc"
echo "Cutadapt is removing poly G reads :)"
date

for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1.fastq.gz)

    if [[ -f "$OUT/${sample}_R2.fastq.gz" ]]; then
        echo "$sample already done, skipping"
        continue
    fi

    echo "Processing: $sample"
    # Minimum length 25 so that reads entirely removed by nextseq-trim are filtered out
    cutadapt --nextseq-trim=20 \
        -j "$THREADS" \
        --minimum-length 25 \
        -o "$OUT/${sample}_R1.fastq.gz" \
        -p "$OUT/${sample}_R2.fastq.gz" \
        "$IN/${sample}_R1.fastq.gz" "$IN/${sample}_R2.fastq.gz"
done

echo "Completed Cutadapt"
date

# FastQC
for fq in "$OUT"/*.fastq.gz; do
    report=$OUT/fastqc/$(basename "$fq" .fastq.gz)_fastqc.html
    [[ -f "$report" ]] && continue
    echo "Running FastQC on $(basename "$fq")..."
    fastqc "$fq" -o "$OUT/fastqc"
done

# MultiQC
conda run -n multiqc_env multiqc -f "$OUT/fastqc" -o "$OUT/fastqc"

echo "All done! Thank you, please come again."
date
