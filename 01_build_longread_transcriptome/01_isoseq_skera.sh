#!/bin/bash

## Step 01 (Iso-Seq): split Kinnex arrays into single transcripts with skera
## Last modified by Linley Sherin 2026-10-22
#
# Runs on the cluster once for all 4 SMRT cells of run SHE27514.
# Each ~19 kb Kinnex read is a concatenation of 8 cDNAs separated by the mas8
# adapters; skera splits them into ~2 kb segmented reads (S-reads).
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/01_isoseq_skera.sh
#
# Reads:  data/01_build_longread_transcriptome/00_hifireads/*.hifi_reads.bcM000N.bam
#           (one per SMRT cell; bcM000N = cell/pool N)
#         data/reference/mas8_primers.fasta
# Writes: data/01_build_longread_transcriptome/01_skera/
#           segmented_cellN.bam      S-reads, input to step 02
#           segmented_cellN.*        skera's other outputs (non_passing, summary)
#           segmented_cellN.log
#
# A cell is skipped if segmented_cellN.bam already exists (e.g. the files
# SickKids provided). Delete a cell's segmented_cellN.* files to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/00_hifireads
OUT=$DATA/01_skera
ADAPTERS=$ROOT/data/reference/mas8_primers.fasta   # Kinnex full-length RNA, 8-mer arrays (MAS-Seq_Adapter_v3)
THREADS=24                                    # speed only, does not change results

[[ -f "$ADAPTERS" ]] || { echo "ERROR: adapter file $ADAPTERS not found" >&2; exit 1; }

shopt -s nullglob
hifi_bams=("$IN"/*.hifi_reads.bcM*.bam)
[[ ${#hifi_bams[@]} -gt 0 ]] || { echo "ERROR: no *.hifi_reads.bcM*.bam files in $IN" >&2; exit 1; }

mkdir -p "$OUT"
echo "Skera is splitting arrays :) [${#hifi_bams[@]} cells]"
date

for bam in "${hifi_bams[@]}"; do
    cell=$(basename "$bam" | sed -E 's/.*\.bcM0*([0-9]+)\.bam$/\1/')
    out_bam=$OUT/segmented_cell$cell.bam

    if [[ -f "$out_bam" ]]; then
        echo "cell$cell already done ($(basename "$out_bam") exists), skipping"
        continue
    fi

    echo "Processing cell$cell: $(basename "$bam")"
    skera split -j "$THREADS" \
        --log-level INFO --log-file "$OUT/segmented_cell$cell.log" \
        "$bam" "$ADAPTERS" "$out_bam"
    echo "Finished cell$cell"
    date
done

echo "Skera is finished :)"
date
