#!/bin/bash

## Step 04 (Iso-Seq): give each FLNC BAM its real sample name (SM tag)
## Last modified by Linley Sherin 2026-10-22
#
# Because we have 24 samples, and Kinnex can only multiplex 12, the same barcodes
# are reused across cells and the SM tag (BioSample1-12) is duplicated. This
# step sets SM to the real sample ID, using updatebam.csv:
#   input_bam,output_bam,new_sample_name
#   03_refine/cell1/flnc_IsoSeqX_bc01_5p--IsoSeqX_3p.bam,04_updatebam/flnc_3FGP.bam,3FGP
# (paths relative to data/01_build_longread_transcriptome/)
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/04_isoseq_updatebam.sh
#
# Reads:  data/01_build_longread_transcriptome/updatebam.csv
#         data/01_build_longread_transcriptome/03_refine/cellN/flnc_*.bam
# Writes: data/01_build_longread_transcriptome/04_updatebam/
#           flnc_<sample>.bam, flnc_<sample>.bam.pbi
#           flnc.fofn                 list of all sample BAMs, input to steps 05 and 06
#
# A sample is skipped if its flnc_<sample>.bam already exists. Delete it to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
CSV=$DATA/updatebam.csv
OUT=$DATA/04_updatebam

[[ -f "$CSV" ]] || { echo "ERROR: $CSV not found" >&2; exit 1; }

mkdir -p "$OUT"
: > "$OUT/flnc.fofn"

echo "Updating BAM sample names..."
date

# Strip Windows line endings, skip the header line
awk 'NR > 1 { sub(/\r$/, ""); print }' "$CSV" |
while IFS=, read -r input_bam output_bam new_sample_name; do
    # Skip empty lines or lines starting with #
    [[ -z "$input_bam" || "$input_bam" =~ ^# ]] && continue

    input=$DATA/$input_bam
    output=$DATA/$output_bam

    # flnc.fofn lists file names only; isoseq finds them next to the fofn
    basename "$output" >> "$OUT/flnc.fofn"

    if [[ -f "$output" ]]; then
        echo "$new_sample_name already done, skipping"
        continue
    fi

    [[ -f "$input" ]] || { echo "ERROR: $input not found. Run step 03 first." >&2; exit 1; }
    echo "Processing: $input_bam -> $output_bam with sample name $new_sample_name"

    # Set SM:<sample> on every @RG line, then replace the header
    samtools view -H "$input" |
        awk -v sm="$new_sample_name" 'BEGIN { FS = OFS = "\t" }
            /^@RG/ { for (i = 1; i <= NF; ++i) if ($i ~ /^SM:/) $i = "SM:" sm }
            { print }' > "$OUT/header.tmp"
    samtools reheader "$OUT/header.tmp" "$input" > "$output"
    rm "$OUT/header.tmp"
    pbindex "$output"
done

echo "Updated BAMs and flnc.fofn written to $OUT"
date
