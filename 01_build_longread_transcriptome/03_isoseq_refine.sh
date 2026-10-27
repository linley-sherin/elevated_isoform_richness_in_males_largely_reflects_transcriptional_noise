#!/bin/bash

## Step 03 (Iso-Seq): remove poly-A tails and concatemers with isoseq refine (FLNC reads)
## Last modified by Linley Sherin 2026-10-22
#
# Runs on the cluster. We don't expect concatemers because skera already split
# the arrays. Refine needs the primer file, so it runs while the reads are still
# in cell-specific folders.
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/03_isoseq_refine.sh
#
# Reads:  data/01_build_longread_transcriptome/02_lima/cellN/output.IsoSeqX_bcNN_5p--IsoSeqX_3p.bam
#         data/reference/IsoSeq_v2_primers_12.fasta
# Writes: data/01_build_longread_transcriptome/03_refine/cellN/
#           flnc_IsoSeqX_bcNN_5p--IsoSeqX_3p.bam     FLNC reads, input to step 04
#           flnc_*.report.csv, flnc_*.filter_summary.report.json
#           flnc_*_refine.log
#
# A barcode is skipped if its flnc_*.bam already exists. Delete its flnc_* files
# to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/02_lima
OUT=$DATA/03_refine
PRIMERS=$ROOT/data/reference/IsoSeq_v2_primers_12.fasta
THREADS=24                                           # speed only, does not change results

[[ -f "$PRIMERS" ]] || { echo "ERROR: primer file $PRIMERS not found" >&2; exit 1; }

shopt -s nullglob
cell_dirs=("$IN"/cell*)
[[ ${#cell_dirs[@]} -gt 0 ]] || { echo "ERROR: no cell folders in $IN. Run step 02 first." >&2; exit 1; }

echo "Refining... [${#cell_dirs[@]} cells]"
date

for cell_dir in "${cell_dirs[@]}"; do
    cell=$(basename "$cell_dir")
    cell_out=$OUT/$cell
    mkdir -p "$cell_out"

    # Only the per-barcode BAMs (skips output.unbarcoded.bam)
    for input in "$cell_dir"/output.IsoSeqX_bc*_5p--IsoSeqX_3p.bam; do
        base_name=$(basename "$input" .bam)
        base_name=${base_name#output.}
        output=$cell_out/flnc_${base_name}.bam

        if [[ -f "$output" ]]; then
            echo "$cell $base_name already done, skipping"
            continue
        fi

        echo "Refining $cell $base_name"
        isoseq refine "$input" "$PRIMERS" "$output" \
            --require-polya \
            -j "$THREADS" \
            --log-level INFO --log-file "$cell_out/flnc_${base_name}_refine.log"
    done

    echo "Refinement for $cell is finished."
    date
done

echo "Refinement is completed :)"
date
