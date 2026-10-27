#!/bin/bash

## Step 04 (short-read transcriptome): quantify short reads against the filtered StringTie transcripts
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally.
# (1) Transcript FASTA (gffread) and t2g map from merged_filtered.gtf
# (2) kb ref + kb count per sample (kb_python wrapper, so the short-read and
#     long-read pipelines match as closely as possible; the abundance matrix is
#     the same as with kallisto quant)
# (3) Mapping rate per sample
#
# Usage (from the project folder):
#   conda activate kallisto_env
#   bash 02_build_shortread_transcriptome/04_kallisto_short.sh
#
# Reads:  data/02_build_shortread_transcriptome/03_filter_stringtie/merged_filtered.gtf
#         data/reference/ncbi_female.fa
#         data/01_build_longread_transcriptome/08_trimmomatic/<sample>_R{1,2}_paired.fastq.gz
# Writes: data/02_build_shortread_transcriptome/
#           data/GTF/05_shortreadpipe.gtf     copy of merged_filtered.gtf, input to step 05
#           04_kallisto/shortreadpipe.fasta, t2g.txt, kallisto_index_short.idx
#           04_kallisto/kallisto_short/<sample>/
#           04_kallisto/maprate_short.txt
#
# (1) is skipped if shortreadpipe.fasta already exists, the index if it already
# exists, and a sample if its run_info.json already exists. Delete them to redo.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
LONG=$ROOT/data/01_build_longread_transcriptome
SHORT=$ROOT/data/02_build_shortread_transcriptome
GTF_CLEAN=$SHORT/03_filter_stringtie/merged_filtered.gtf
GENOME=$ROOT/data/reference/ncbi_female.fa
READS=$LONG/08_trimmomatic
OUT=$SHORT/04_kallisto
FASTA=$OUT/shortreadpipe.fasta
T2G=$OUT/t2g.txt
INDEX=$OUT/kallisto_index_short.idx
GFFREAD=~/anaconda3/bin/gffread    # in conda base, not in kallisto_env

mkdir -p "$OUT/kallisto_short" "$ROOT/data/GTF"

## (1) Transcript FASTA and t2g map
if [[ -f "$FASTA" ]]; then
    echo "shortreadpipe.fasta already exists, skipping FASTA and t2g"
else
    [[ -f "$GTF_CLEAN" ]] || { echo "ERROR: $GTF_CLEAN not found. Run step 03 first." >&2; exit 1; }

    echo "Generating spliced transcript FASTA..."
    "$GFFREAD" "$GTF_CLEAN" -g "$GENOME" -w "$FASTA"

    echo "Generating t2g mapping..."
    # Split the attributes column by double quotes and pick out transcript_id and gene_id
    awk -F'\t' '$3=="transcript" {
        split($9, a, "\"");
        tid = ""; gid = "";
        for (i = 1; i <= length(a); i++) {
            if (a[i] ~ /transcript_id/) tid = a[i+1];
            if (a[i] ~ /gene_id/) gid = a[i+1];
        }
        if (tid != "" && gid != "") print tid "\t" gid;
    }' "$GTF_CLEAN" > "$T2G"

    # Final GTF copy for your records
    cp "$GTF_CLEAN" "$ROOT/data/GTF/05_shortreadpipe.gtf"
fi

## (2) kallisto index and quantification
if [[ -f "$INDEX" ]]; then
    echo "Kallisto index already exists, skipping index creation"
else
    echo "Creating kallisto index..."
    kb ref "$FASTA" \
        -i "$INDEX" \
        --make-unique \
        --workflow=custom \
        --tmp "$OUT/tmp_short"
fi

shopt -s nullglob
r1_files=("$READS"/*_R1_paired.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1_paired.fastq.gz files in $READS. Run long-read step 08 first." >&2; exit 1; }

for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1_paired.fastq.gz)
    sample_out=$OUT/kallisto_short/$sample

    if [[ -f "$sample_out/run_info.json" ]]; then
        echo "$sample already done, skipping"
        continue
    fi

    echo "Running kallisto quantification for $sample..."
    mkdir -p "$sample_out"
    kb count -x bulk \
        -i "$INDEX" \
        -o "$sample_out" \
        -g "$T2G" \
        --parity=paired \
        --strand=unstranded \
        --tcc \
        --matrix-to-directories \
        "$READS/${sample}_R1_paired.fastq.gz" \
        "$READS/${sample}_R2_paired.fastq.gz" \
        --bootstraps 100 \
        --tmp "$OUT/tmp_short"
done

## (3) Mapping rate per sample
echo "Finding mapping rate for short reads..."
: > "$OUT/maprate_short.txt"
for run_info in "$OUT"/kallisto_short/*/run_info.json; do
    sample=$(basename "$(dirname "$run_info")")
    processed_reads=$(jq -r '.n_processed' "$run_info")
    mapped_reads=$(jq -r '.n_pseudoaligned' "$run_info")
    pseudoaligned=$(jq -r '.p_pseudoaligned' "$run_info")
    echo "$sample processed_reads: $processed_reads mapped_reads: $mapped_reads p_pseudoaligned: $pseudoaligned" >> "$OUT/maprate_short.txt"
done
echo "Mapping rates saved in $OUT/maprate_short.txt"

echo "Kallisto is finished :)"
date
