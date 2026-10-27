#!/bin/bash

## Step 15: quantify long reads (FLNC) against the final non-redundant transcriptome
## Last modified by Linley Sherin 2026-10-23
#
# Runs locally.
# (1) t2g map (transcript -> gene) from the non-redundant transcript FASTA
# (2) FLNC BAMs -> FASTQ (only for samples without a FASTQ yet)
# (3) kb ref + kb count --long per sample (kb_python wrapper, so the short-read and
#     long-read pipelines match as closely as possible). Needs the kallisto build for
#     k-mer size 64 (kallisto_LongKmer_NoOpt-v0.51.1), in tools/kallisto/ (project folder).
# (4) Mapping rate per sample
#
# Usage (from the project folder):
#   conda activate kallisto_env
#   bash 01_build_longread_transcriptome/15_kallisto_long.sh
#
# Reads:  data/01_build_longread_transcriptome/14_remove_redundancy/nonredundant_transcripts.fasta
#         data/01_build_longread_transcriptome/04_updatebam/flnc_<sample>.bam (or .fastq)
#         tools/kallisto/kallisto_mac-v0.51.1_kmer64
# Writes: data/01_build_longread_transcriptome/04_updatebam/flnc_<sample>.fastq (if missing)
#         data/01_build_longread_transcriptome/15_kallisto/
#           t2g.txt, kallisto_index_long.idx
#           kallisto_long/flnc_<sample>.fastq/   kb output per sample, input to step 16
#           maprate_long.txt
#
# The t2g and index are skipped if they already exist, and a sample if its
# run_info.json already exists. Delete them to redo.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
FASTA=$DATA/14_remove_redundancy/nonredundant_transcripts.fasta
FLNC=$DATA/04_updatebam
KALLISTO=$ROOT/tools/kallisto/kallisto_mac-v0.51.1_kmer64
OUT=$DATA/15_kallisto
T2G=$OUT/t2g.txt
INDEX=$OUT/kallisto_index_long.idx
THREADS=40          # speed only, does not change results

[[ -f "$FASTA" ]] || { echo "ERROR: $FASTA not found. Run step 14 first." >&2; exit 1; }
[[ -x "$KALLISTO" ]] || { echo "ERROR: kallisto k64 build not found at $KALLISTO" >&2; exit 1; }

mkdir -p "$OUT/kallisto_long"

## (1) t2g map: gene name = transcript name without the version number (.1, .2, ...)
if [[ -f "$T2G" ]]; then
    echo "t2g.txt already exists, skipping"
else
    echo "Creating t2g map..."
    grep '^>' "$FASTA" | sed 's/^>//' | awk '{ gene = $1; sub(/\.[0-9]+$/, "", gene); print $1 "\t" gene }' > "$T2G"
fi

## (2) FLNC BAMs -> FASTQ for kallisto
shopt -s nullglob
for bam in "$FLNC"/flnc_*.bam; do
    fastq=${bam%.bam}.fastq
    [[ -f "$fastq" ]] && continue
    echo "Converting $(basename "$bam") to FASTQ..."
    samtools fastq "$bam" > "$fastq"
done

fastqs=("$FLNC"/flnc_*.fastq)
[[ ${#fastqs[@]} -gt 0 ]] || { echo "ERROR: no flnc_*.fastq files in $FLNC. Run step 04 first." >&2; exit 1; }

## (3) kallisto index and quantification
if [[ -f "$INDEX" ]]; then
    echo "Kallisto index already exists, skipping index creation"
else
    echo "Creating kallisto index..."
    kb ref "$FASTA" \
        -i "$INDEX" \
        --make-unique \
        --workflow=custom \
        -k 63 \
        --kallisto "$KALLISTO" \
        --tmp "$OUT/tmp_long"
fi

for fastq in "${fastqs[@]}"; do
    sample_out=$OUT/kallisto_long/$(basename "$fastq")

    if [[ -f "$sample_out/run_info.json" ]]; then
        echo "$(basename "$fastq") already done, skipping"
        continue
    fi

    # 100 bootstraps so uncertainty can be calculated later
    echo "Running kallisto quantification for $(basename "$fastq")..."
    mkdir -p "$sample_out"
    kb count -x bulk \
        -i "$INDEX" \
        -o "$sample_out" \
        -g "$T2G" \
        --parity=single \
        --strand=unstranded \
        --tcc \
        --matrix-to-directories \
        --long \
        --platform PacBio \
        "$fastq" \
        --kallisto "$KALLISTO" \
        --bootstraps 100 \
        -t "$THREADS" \
        --tmp "$OUT/tmp_long"
done

## (4) Mapping rate per sample
echo "Finding mapping rate for long reads..."
: > "$OUT/maprate_long.txt"
for run_info in "$OUT"/kallisto_long/flnc_*/run_info.json; do
    sample=$(basename "$(dirname "$run_info")")
    sample=${sample#flnc_}
    sample=${sample%.fastq}
    processed_reads=$(jq -r '.n_processed' "$run_info")
    mapped_reads=$(jq -r '.n_pseudoaligned' "$run_info")
    pseudoaligned=$(jq -r '.p_pseudoaligned' "$run_info")
    echo "$sample processed_reads: $processed_reads mapped_reads: $mapped_reads p_pseudoaligned: $pseudoaligned" >> "$OUT/maprate_long.txt"
done
echo "Mapping rates saved in $OUT/maprate_long.txt"

echo "Kallisto is finished :)"
date
