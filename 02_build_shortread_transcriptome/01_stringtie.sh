#!/bin/bash

## Step 01 (short-read transcriptome): assemble transcripts per sample with StringTie, then merge
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally. Uses the STAR BAMs from the long-read section (step 09, NCBI genome),
# guided by the NCBI annotation (ncbi_female_stringtie.gtf, edited to work with StringTie).
#
# Usage (from the project folder; StringTie is called directly, no env):
#   bash 02_build_shortread_transcriptome/01_stringtie.sh
#
# Reads:  data/01_build_longread_transcriptome/09_star/<sample>_Aligned.sortedByCoord.out.bam
#         data/reference/ncbi_female_stringtie.gtf
# Writes: data/02_build_shortread_transcriptome/01_stringtie/
#           per_sample/<sample>.gtf    one assembly per sample
#           merge/gtf_list.txt
#           merge/merged.gtf           all samples merged, input to steps 02 and 03
#
# A sample is skipped if its per_sample/<sample>.gtf already exists, and the merge
# if merged.gtf already exists. Delete them to redo either part.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
LONG=$ROOT/data/01_build_longread_transcriptome
BAM_DIR=$LONG/09_star
REF_GTF=$ROOT/data/reference/ncbi_female_stringtie.gtf
OUT=$ROOT/data/02_build_shortread_transcriptome/01_stringtie
THREADS=20          # speed only, does not change results

[[ -f "$REF_GTF" ]] || { echo "ERROR: $REF_GTF not found" >&2; exit 1; }

shopt -s nullglob
bams=("$BAM_DIR"/*_Aligned.sortedByCoord.out.bam)
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: no STAR BAMs in $BAM_DIR. Run long-read step 09 first." >&2; exit 1; }

mkdir -p "$OUT/per_sample" "$OUT/merge"

# Initial assembly: one transcriptome per sample
for bam in "${bams[@]}"; do
    sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
    if [[ -f "$OUT/per_sample/$sample.gtf" ]]; then
        echo "$sample already done, skipping"
        continue
    fi
    echo "Processing sample: $sample"
    stringtie "$bam" \
        -G "$REF_GTF" \
        -o "$OUT/per_sample/$sample.gtf" \
        -p "$THREADS"
done

# Merge assemblies into one merged.gtf
if [[ -f "$OUT/merge/merged.gtf" ]]; then
    echo "merged.gtf already exists, skipping merge"
else
    echo "Merging all GTFs..."
    ls "$OUT/per_sample/"*.gtf > "$OUT/merge/gtf_list.txt"
    stringtie --merge \
        -G "$REF_GTF" \
        -o "$OUT/merge/merged.gtf" \
        -p "$THREADS" \
        "$OUT/merge/gtf_list.txt"
fi

echo "Done! Results are in $OUT"
date
