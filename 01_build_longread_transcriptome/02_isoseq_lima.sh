#!/bin/bash

## Step 02 (Iso-Seq): remove primers and demultiplex S-reads by sample barcode with lima
## Last modified by Linley Sherin 2026-10-22
#
# Runs on the cluster. Each cell is done separately because the same
# barcodes (bc01-12) are reused across cells.
#
# NOTE: to get the original primer file from pacbio:
# wget https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/02_isoseq_lima.sh
#
# Reads:  data/01_build_longread_transcriptome/01_skera/segmented_cellN.bam
#         data/reference/IsoSeq_v2_primers_12.fasta
# Writes: data/01_build_longread_transcriptome/02_lima/cellN/
#           output.IsoSeqX_bcNN_5p--IsoSeqX_3p.bam   one BAM per barcode
#           output.lima.*                            reports (unbarcoded reads are not kept)
#           cellN_lima.log
#
# A cell is skipped if output.lima.summary already exists (lima writes it
# last). Delete a cell's folder to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/01_skera
OUT=$DATA/02_lima
PRIMERS=$ROOT/data/reference/IsoSeq_v2_primers_12.fasta   # IsoSeqX bc01-12 5' primers + 3' primer
THREADS=24                                           # speed only, does not change results

[[ -f "$PRIMERS" ]] || { echo "ERROR: primer file $PRIMERS not found" >&2; exit 1; }

shopt -s nullglob
# Keep only segmented_cellN.bam, not skera's other outputs (e.g. segmented_cellN.non_passing.bam)
segmented_bams=()
for bam in "$IN"/segmented_cell*.bam; do
    [[ $(basename "$bam") =~ ^segmented_cell[0-9]+\.bam$ ]] && segmented_bams+=("$bam")
done
[[ ${#segmented_bams[@]} -gt 0 ]] || { echo "ERROR: no segmented_cell*.bam files in $IN. Run step 01 first." >&2; exit 1; }

echo "Starting lima! [${#segmented_bams[@]} cells]"
date

for bam in "${segmented_bams[@]}"; do
    cell=$(basename "$bam" .bam | sed 's/^segmented_//')   # e.g. cell1
    cell_out=$OUT/$cell

    if [[ -f "$cell_out/output.lima.summary" ]]; then
        echo "$cell already done, skipping"
        continue
    fi

    echo "Processing $cell"
    mkdir -p "$cell_out"
    lima "$bam" "$PRIMERS" "$cell_out/output.bam" \
        --isoseq \
        --peek-guess \
        --overwrite-biosample-names \
        -j "$THREADS" \
        --log-level INFO --log-file "$cell_out/${cell}_lima.log"
    echo "Processing for $cell is finished."
    date
done

echo "Lima is finished :)"
date
