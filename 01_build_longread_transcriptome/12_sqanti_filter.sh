#!/bin/bash

## Step 12: SQANTI3 rules filter
## Last modified by Linley Sherin 2026-10-23
#
# Runs locally. Rules are in filter_new_filterstruc.json:
#   - Keep all FSM isoform models without evidence of intrapriming (%A downstream of TTS >60,
#     ie. 12/20 bp; indicates 3' end is potentially unreliable)
#   - Keep non-FSM isoform models if they meet the below criteria:
#     (1) no evidence of intrapriming (see above)
#     (2) no junction labeled RT-switching (RTS)
#     (3) all junctions canonical OR have short-read cov >10 (-c threshold)
#
# Usage (from the project folder):
#   conda activate SQANTI3.env
#   bash 01_build_longread_transcriptome/12_sqanti_filter.sh
#
# Reads:  data/01_build_longread_transcriptome/11_sqanti_qc/output/SQANTI_output_classification.txt,
#           SQANTI_output_corrected.fasta, SQANTI_output_corrected.gtf
#         data/01_build_longread_transcriptome/12_sqanti_filter/filter_new_filterstruc.json
# Writes: data/01_build_longread_transcriptome/12_sqanti_filter/
#           SQANTI_filtered.filtered.gtf, .fasta    input to steps 13 and 14
#           SQANTI_filtered_RulesFilter_result_classification.txt, ..._filtering_reasons.txt, report
#         data/GTF/02_SQANTI_filtered.gtf           copy of SQANTI_filtered.filtered.gtf
#
# Skipped if SQANTI_filtered.filtered.gtf already exists. Delete it to redo this step.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/data/01_build_longread_transcriptome
QC=$DATA/11_sqanti_qc/output
OUT=$DATA/12_sqanti_filter
RULES=$OUT/filter_new_filterstruc.json
SQANTI_FILTER=~/Documents/programs/SQANTI3-5.2.2/sqanti3_filter.py

for f in "$QC/SQANTI_output_classification.txt" "$QC/SQANTI_output_corrected.fasta" "$QC/SQANTI_output_corrected.gtf" "$RULES"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done

if [[ -f "$OUT/SQANTI_filtered.filtered.gtf" ]]; then
    echo "SQANTI_filtered.filtered.gtf already exists, skipping"
    exit 0
fi

echo "Running SQANTI3 rules filter..."
date

python "$SQANTI_FILTER" rules "$QC/SQANTI_output_classification.txt" `# classification file from sqanti3` \
    --isoforms "$QC/SQANTI_output_corrected.fasta" `# fasta file with isoform models to filter` \
    --gtf "$QC/SQANTI_output_corrected.gtf" `# gtf file with isoform models to filter` \
    -o SQANTI_filtered \
    -d "$OUT" \
    -j "$RULES" `# json with rules`

cp "$OUT/SQANTI_filtered.filtered.gtf" "$ROOT/data/GTF/02_SQANTI_filtered.gtf"

echo "SQANTI3 filter is finished :)"
date
