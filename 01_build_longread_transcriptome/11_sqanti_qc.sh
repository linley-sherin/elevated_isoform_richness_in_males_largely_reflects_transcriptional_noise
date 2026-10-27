#!/bin/bash

## Step 11: SQANTI3 QC of the collapsed isoforms, with short-read support
## Last modified by Linley Sherin 2026-10-22
#
# STAR (step 09) and kallisto (step 10) were already run, so their output is
# given to SQANTI3 instead of having SQANTI3 rerun them (no --short_reads):
#   -c        STAR splice junctions (09_star/*_SJ.out.tab)
#   --SR_bam  STAR BAMs (listed in STAR_bams.fofn, made below)
#   --expression  kallisto tables, copied to kallisto_abund/<sample>.tsv below
#             with the IDs trimmed to match the GTF (PB.1.1|PB.1.1:...|transcript/... -> PB.1.1)
#
# NOTE: Can only run on mank03 and mank05
# Must initiate screen with csh (not bash):
#   screen -S [session_name] csh
#
# Usage (on mank03/mank05, from the project folder):
#   bash 01_build_longread_transcriptome/11_sqanti_qc.sh
#
# Reads:  data/01_build_longread_transcriptome/06_mapcollapse/collapsed.gff, collapsed.flnc_count.txt
#         data/01_build_longread_transcriptome/09_star/<sample>_SJ.out.tab, <sample>_Aligned.sortedByCoord.out.bam
#         data/01_build_longread_transcriptome/10_kallisto/<sample>/abundance.tsv
#         data/reference/ncbi_female.fa, ncbi_female_sorted.gtf
# Writes: data/01_build_longread_transcriptome/11_sqanti_qc/
#           kallisto_abund/<sample>.tsv   kallisto tables for SQANTI3
#           STAR_bams.fofn                list of STAR BAMs for SQANTI3
#           output/SQANTI_output_*        classification, junctions, corrected GTF, report
#           output/sqanti_terminal.log
#         data/GTF/01_SQANTI_output.gtf  copy of output/SQANTI_output_corrected.gtf
#
# Skipped if output/SQANTI_output_classification.txt already exists. Delete
# output/ to redo it.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
ISOFORMS=$DATA/06_mapcollapse/collapsed.gff
FL_COUNTS=$DATA/06_mapcollapse/collapsed.flnc_count.txt
STAR_DIR=$DATA/09_star
KALLISTO_DIR=$DATA/10_kallisto
ANNOTATION=$ROOT/data/reference/ncbi_female_sorted.gtf
GENOME=$ROOT/data/reference/ncbi_female.fa
OUT=$DATA/11_sqanti_qc
CPUS=48             # speed only, does not change results
CHUNKS=12           # speed only, does not change results

for f in "$ISOFORMS" "$FL_COUNTS" "$ANNOTATION" "$GENOME"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done

if [[ -f "$OUT/output/SQANTI_output_classification.txt" ]]; then
    echo "SQANTI3 output already exists, skipping"
    exit 0
fi

mkdir -p "$OUT/output" "$OUT/kallisto_abund"

# STAR BAMs -> fofn
shopt -s nullglob
star_bams=("$STAR_DIR"/*_Aligned.sortedByCoord.out.bam)
[[ ${#star_bams[@]} -gt 0 ]] || { echo "ERROR: no STAR BAMs in $STAR_DIR. Run step 09 first." >&2; exit 1; }
printf "%s\n" "${star_bams[@]}" > "$OUT/STAR_bams.fofn"

# kallisto tables -> kallisto_abund/<sample>.tsv, IDs trimmed at the first "|"
kallisto_tables=("$KALLISTO_DIR"/*/abundance.tsv)
[[ ${#kallisto_tables[@]} -gt 0 ]] || { echo "ERROR: no abundance.tsv files in $KALLISTO_DIR/*/. Run step 10 first." >&2; exit 1; }
for tsv in "${kallisto_tables[@]}"; do
    sample=$(basename "$(dirname "$tsv")")
    awk 'BEGIN { FS = OFS = "\t" } NR > 1 { sub(/\|.*/, "", $1) } { print }' "$tsv" > "$OUT/kallisto_abund/$sample.tsv"
done

echo "Running SQANTI3 QC [${#star_bams[@]} STAR BAMs, ${#kallisto_tables[@]} kallisto tables]..."
date

apptainer run /Linux/containers/sqanti3.sif sqanti3_qc.py \
    "$ISOFORMS" `# long-read transcriptome` \
    "$ANNOTATION" `# female guppy annotation` \
    "$GENOME" `# female guppy genome` \
    -o SQANTI_output `# output prefix` \
    -d "$OUT/output" `# output directory` \
    -fl "$FL_COUNTS" `# long-read counts (flnc)` \
    -c "$STAR_DIR/*_SJ.out.tab" `# splice junction coverage from STAR` \
    --SR_bam "$OUT/STAR_bams.fofn" `# mapped short-read BAMs from STAR` \
    --expression "$OUT/kallisto_abund/" `# directory of abundance tables from kallisto` \
    --saturation `# include saturation curves in report` \
    --skipORF \
    --cpus "$CPUS" \
    --chunks "$CHUNKS" \
    --report both > "$OUT/output/sqanti_terminal.log" 2>&1

cp "$OUT/output/SQANTI_output_corrected.gtf" "$ROOT/data/GTF/01_SQANTI_output.gtf"

echo "SQANTI3 QC is finished :)"
date
