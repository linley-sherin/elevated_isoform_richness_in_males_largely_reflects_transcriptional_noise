#!/bin/bash

## Step 13: ORF prediction with TD2 on the SQANTI-filtered isoforms, then best ORF per transcript
## Last modified by Linley Sherin 2026-10-23
#
# Runs locally.
# (1) Transcript FASTA + GFF3 from the SQANTI-filtered GTF (TransDecoder utils)
# (2) TD2.LongOrfs + TD2.Predict (-P 0.9)
# (3) ORFs mapped back to genome coordinates (transcripts_with_cds.gff3)
# (4) Best ORF per transcript (lib/parse_td2.R)
#
# Usage (from the project folder):
#   conda activate transdecoder_env
#   bash 01_build_longread_transcriptome/13_transdecoder.sh
#
# Reads:  data/GTF/02_SQANTI_filtered.gtf
#         data/reference/ncbi_female.fa
#         tools/TransDecoder-v5.7.1/util/
# Writes: data/01_build_longread_transcriptome/13_transdecoder/
#           transcripts.fasta, transcripts.gff3
#           transcripts.fasta.TD2.*, longest_orfs.*
#           transcripts_with_cds.gff3
#           parsed_results.csv          best ORF per transcript, input to step 14
#           transdecoder2_<date>.log
#
# TD2 (1-3) is skipped if transcripts_with_cds.gff3 already exists, and (4) if
# parsed_results.csv already exists. Delete them to redo either part.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
GTF=$ROOT/data/GTF/02_SQANTI_filtered.gtf
GENOME=$ROOT/data/reference/ncbi_female.fa
TD_UTIL=$ROOT/tools/TransDecoder-v5.7.1/util
OUT=$ROOT/data/01_build_longread_transcriptome/13_transdecoder
RSCRIPT=/usr/local/bin/Rscript      # system R (transdecoder_env has its own R without the packages)

export PATH=/Library/Frameworks/Python.framework/Versions/3.11/bin/TD2/util:$PATH

[[ -f "$GTF" ]] || { echo "ERROR: $GTF not found. Run step 12 first." >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: genome $GENOME not found" >&2; exit 1; }

mkdir -p "$OUT"

if [[ -f "$OUT/transcripts_with_cds.gff3" ]]; then
    echo "transcripts_with_cds.gff3 already exists, skipping TD2"
else
    # === Setup logging ===
    LOGFILE="$OUT/transdecoder2_$(date +%Y%m%d_%H%M%S).log"
    exec > >(tee -i "$LOGFILE") 2>&1
    echo "Starting TransDecoder analysis at $(date)"

    # (1) Transcript FASTA directly from genomic coordinates, and GTF -> GFF3
    "$TD_UTIL/gtf_genome_to_cdna_fasta.pl" "$GTF" "$GENOME" > "$OUT/transcripts.fasta"
    "$TD_UTIL/gtf_to_alignment_gff3.pl" "$GTF" > "$OUT/transcripts.gff3"

    # (2) TD2 runs inside the output folder
    cd "$OUT"
    if [[ -d transcripts ]]; then
        echo "TD2.LongOrfs output folder 'transcripts' already exists, skipping"
    else
        echo "Running TD2.LongOrfs (relaxed)..."
        TD2.LongOrfs -t transcripts.fasta -v -S
    fi

    echo "Running TD2.Predict..."
    TD2.Predict -t transcripts.fasta -v -O ./transcripts -P 0.9

    echo "Moving final output files up one level..."
    mv transcripts/* ./
    rmdir transcripts
    cd "$ROOT"

    # (3) Map ORFs back to the genome (fixes strand mismatches)
    "$TD_UTIL/cdna_alignment_orf_to_genome_orf.pl" \
        "$OUT/transcripts.fasta.TD2.gff3" \
        "$OUT/transcripts.gff3" \
        "$OUT/transcripts.fasta" > "$OUT/transcripts_with_cds.gff3"

    echo "TransDecoder finished at $(date)"
fi

# (4) Best ORF per transcript
if [[ -f "$OUT/parsed_results.csv" ]]; then
    echo "parsed_results.csv already exists, skipping parse"
else
    echo "Parsing TD2 results..."
    "$RSCRIPT" "$ROOT/01_build_longread_transcriptome/lib/parse_td2.R"
fi

echo "Step 13 is finished :)"
date
