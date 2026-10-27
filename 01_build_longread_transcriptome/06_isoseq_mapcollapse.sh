#!/bin/bash

## Step 06 (Iso-Seq): map clusters to the genome (pbmm2) and collapse into unique isoforms
## Last modified by Linley Sherin 2026-10-22
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/06_isoseq_mapcollapse.sh
#
# Reads:  data/01_build_longread_transcriptome/05_cluster/clustered.bam
#         data/01_build_longread_transcriptome/04_updatebam/flnc.fofn (and the BAMs it lists)
#         data/reference/ncbi_female.fa
# Writes: data/01_build_longread_transcriptome/06_mapcollapse/
#           mapped.bam, mapped.log
#           collapsed.gff, collapsed.fasta, collapsed.abundance.txt, ...
#           collapse.log
#
# Mapping is skipped if mapped.bam already exists, and collapsing if
# collapsed.gff already exists. Delete them to redo either part.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
GENOME=$ROOT/data/reference/ncbi_female.fa
CLUSTERED=$DATA/05_cluster/clustered.bam
FOFN=$DATA/04_updatebam/flnc.fofn
OUT=$DATA/06_mapcollapse
PBMM2_THREADS=40      # speed only, does not change results
COLLAPSE_THREADS=48   # speed only, does not change results

[[ -f "$GENOME" ]] || { echo "ERROR: genome $GENOME not found" >&2; exit 1; }
[[ -f "$CLUSTERED" ]] || { echo "ERROR: $CLUSTERED not found. Run step 05 first." >&2; exit 1; }
[[ -f "$FOFN" ]] || { echo "ERROR: $FOFN not found. Run step 04 first." >&2; exit 1; }

mkdir -p "$OUT"

# Align clusters to the genome
if [[ -f "$OUT/mapped.bam" ]]; then
    echo "mapped.bam already exists, skipping mapping"
else
    echo "Mapping clustered transcripts to $(basename "$GENOME")..."
    date
    pbmm2 align "$GENOME" "$CLUSTERED" "$OUT/mapped.bam" \
        --preset ISOSEQ \
        --sort \
        -j "$PBMM2_THREADS" \
        --log-level INFO --log-file "$OUT/mapped.log"
fi

# Collapse mapped transcripts into unique isoforms
if [[ -f "$OUT/collapsed.gff" ]]; then
    echo "collapsed.gff already exists, skipping collapse"
else
    echo "Collapsing isoforms..."
    date
    isoseq collapse "$OUT/mapped.bam" "$FOFN" "$OUT/collapsed.gff" \
        --max-fuzzy-junction 0 \
        --do-not-collapse-extra-5exons \
        -j "$COLLAPSE_THREADS" \
        --log-level INFO --log-file "$OUT/collapse.log"
fi

echo "Map + collapse is finished :)"
date
