#!/bin/bash

## Step 09 (Illumina): STAR index + mapping, for splice junctions and coverage in SQANTI3
## Last modified by Linley Sherin 2026-10-22
#
# Runs locally. SQANTI3 can run STAR itself, but I had issues with it, so this
# runs STAR separately with the same options SQANTI3 uses.
#
# The index is built without an annotation GTF (the STAR manual encourages one,
# but SQANTI3 doesn't want the splice junctions biased towards any annotation).
# It uses the NCBI genome, the same one as Iso-Seq mapping (step 06) and SQANTI3,
# so the chromosome names match (NC_024331.1, ...).
#
# Usage (from the project folder; STAR is called directly, no env):
#   bash 01_build_longread_transcriptome/09_illumina_star.sh
#
# Reads:  data/01_build_longread_transcriptome/08_trimmomatic/<sample>_R{1,2}_paired.fastq.gz
#         data/reference/ncbi_female.fa
# Writes: data/01_build_longread_transcriptome/09_star/
#           index/                                genome index
#           <sample>_Aligned.sortedByCoord.out.bam
#           <sample>_SJ.out.tab                   splice junctions
#           <sample>_Log.final.out, ...
#
# The index is skipped if index/SA already exists, and a sample is skipped
# if its <sample>_SJ.out.tab already exists. Delete them to redo either part.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
IN=$DATA/08_trimmomatic
GENOME=$ROOT/data/reference/ncbi_female.fa
OUT=$DATA/09_star
INDEX=$OUT/index
INDEX_THREADS=4     # speed only, does not change results
MAP_THREADS=40      # speed only, does not change results

[[ -f "$GENOME" ]] || { echo "ERROR: genome $GENOME not found" >&2; exit 1; }

shopt -s nullglob
r1_files=("$IN"/*_R1_paired.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: no *_R1_paired.fastq.gz files in $IN. Run step 08 first." >&2; exit 1; }

mkdir -p "$INDEX"

## Create a STAR index of the guppy genome
if [[ -f "$INDEX/SA" ]]; then
    echo "STAR index already exists, skipping index creation"
else
    echo "Building STAR index from $(basename "$GENOME")..."
    date
    STAR --runMode genomeGenerate \
        --genomeDir "$INDEX" \
        --genomeFastaFiles "$GENOME" \
        --outTmpDir "$INDEX/STAR_tmp" \
        --genomeSAindexNbases 13 \
        --runThreadN "$INDEX_THREADS"
fi

## Map each sample to the genome
for r1 in "${r1_files[@]}"; do
    sample=$(basename "$r1" _R1_paired.fastq.gz)
    prefix=$OUT/${sample}_

    if [[ -f "${prefix}SJ.out.tab" ]]; then
        echo "$sample already done, skipping"
        continue
    fi

    echo "Mapping: $sample"
    date
    STAR --genomeDir "$INDEX" \
        --readFilesIn "$IN/${sample}_R1_paired.fastq.gz" "$IN/${sample}_R2_paired.fastq.gz" \
        --readFilesCommand gzcat \
        --outFileNamePrefix "$prefix" \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterType BySJout \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 1 \
        --genomeLoad NoSharedMemory \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --runThreadN "$MAP_THREADS"
    echo "Finished mapping $sample"
done

echo "Completed STAR mapping."
date
