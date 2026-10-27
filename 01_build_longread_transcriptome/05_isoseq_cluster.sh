#!/bin/bash

## Step 05 (Iso-Seq): cluster FLNC reads from all samples together with isoseq cluster2
## Last modified by Linley Sherin 2026-10-22
#
# Usage (on the cluster, from the project folder):
#   conda activate pacbio
#   bash 01_build_longread_transcriptome/05_isoseq_cluster.sh
#
# Reads:  data/01_build_longread_transcriptome/04_updatebam/flnc.fofn (and the BAMs it lists)
# Writes: data/01_build_longread_transcriptome/05_cluster/
#           clustered.bam             input to step 06
#           clustered.cluster_report.csv
#           cluster.log
#
# Skipped if clustered.bam already exists. Delete it to redo this step.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
FOFN=$DATA/04_updatebam/flnc.fofn
OUT=$DATA/05_cluster
THREADS=48          # speed only, does not change results
SORT_THREADS=24     # speed only, does not change results

[[ -f "$FOFN" ]] || { echo "ERROR: $FOFN not found. Run step 04 first." >&2; exit 1; }

if [[ -f "$OUT/clustered.bam" ]]; then
    echo "$OUT/clustered.bam already exists, skipping"
    exit 0
fi

echo "Clustering $(wc -l < "$FOFN" | tr -d ' ') samples..."
date

mkdir -p "$OUT"
isoseq cluster2 "$FOFN" "$OUT/clustered.bam" \
    -j "$THREADS" \
    --sort-threads "$SORT_THREADS" \
    --log-level INFO --log-file "$OUT/cluster.log"

echo "Clustering is finished :)"
date
