#!/bin/bash

## Step 10 (Illumina): kallisto quantification against the Iso-Seq transcripts (for SQANTI3)
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally.
#
# Usage (from the project folder):
#   conda activate kallisto_env
#   bash 01_build_longread_transcriptome/10_illumina_kallisto.sh
#
# Reads:  data/01_build_longread_transcriptome/08_trimmomatic/<sample>_R{1,2}_paired.fastq.gz
#         data/01_build_longread_transcriptome/06_mapcollapse/collapsed.fasta
# Writes: data/01_build_longread_transcriptome/10_kallisto/
#           kallisto_collapsed.idx
#           <sample>/abundance.tsv, abundance.h5, run_info.json
#
# The index is skipped if it already exists, and a sample is skipped if its
# abundance.h5 already exists. Delete them to redo either part.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/08_trimmomatic
TRANSCRIPTS=$DATA/06_mapcollapse/collapsed.fasta
OUT=$DATA/10_kallisto
INDEX=$OUT/kallisto_collapsed.idx
THREADS=40          # speed only, does not change results

shopt -s nullglob
r1_files=("$IN"/*_R1_paired.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1_paired.fastq.gz files in $IN. Run step 08 first." >&2; exit 1; }

mkdir -p "$OUT"

# Create the kallisto index
if [[ -f "$INDEX" ]]; then
    echo "Kallisto index already exists, skipping index creation"
else
    [[ -f "$TRANSCRIPTS" ]] || { echo "ERROR: $TRANSCRIPTS not found. Run step 06 first." >&2; exit 1; }
    echo "Creating kallisto index..."
    kallisto index -i "$INDEX" "$TRANSCRIPTS" --make-unique
fi

# Quantify each sample
for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1_paired.fastq.gz)
    sample_out=$OUT/$sample

    if [[ -f "$sample_out/abundance.h5" ]]; then
        echo "$sample already done, skipping"
        continue
    fi

    echo "Running kallisto quantification for $sample..."
    mkdir -p "$sample_out"
    kallisto quant -i "$INDEX" \
        -o "$sample_out" \
        -b 100 \
        -t "$THREADS" \
        "$IN/${sample}_R1_paired.fastq.gz" "$IN/${sample}_R2_paired.fastq.gz"
done

echo "Kallisto is finished :)"
date
