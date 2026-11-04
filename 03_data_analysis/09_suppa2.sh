#!/bin/bash

## Step 09: SUPPA2 alternative splicing events and differential splicing (female vs male)
## Last modified by Linley Sherin 2026-11-03
#
# Runs locally. For each tissue (gonad, tail) and each dataset (short, long reads):
# (1) Reformat the TMM-TPM matrix for SUPPA2
# (2) Generate events from the GTF and calculate PSI per event
# (3) Split TPM and PSI into females and males and run diffSplice
# (4) Gather all events into one summary per tissue, then combine tissues
#
# Usage (from the project folder; SUPPA2 is in tools/):
#   bash 03_data_analysis/09_suppa2.sh
#
# Reads:  data/02_build_shortread_transcriptome/05_edgeR/0[12]_*_TMMTPMmatrix_short.csv
#         data/01_build_longread_transcriptome/16_edgeR_filterexp/0[12]_*_TMMTPMmatrix.csv
#         data/GTF/06_shortreadpipe_filtered_by_exp_short.gtf
#         data/GTF/04_filtered_by_exp_long.gtf
# Writes: data/03_data_analysis/09_suppa2/
#           {short,long}_tpm_{gonad,tail}.tsv    TPM in SUPPA2 format
#           {short,long}_{gonad,tail}/           events, psi, diff (dpsi) per dataset and tissue
#           events_summary_{gonad,tail}.tsv, events_summary_combined.tsv
#
# A tissue is skipped if its events_summary_<tissue>.tsv already exists.
# Delete it (and the tissue's folders) to redo that tissue.

set -euo pipefail

# Paths are relative to the project root (one folder up from this script)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

SUPPA_PATH="tools/SUPPA/suppa.py"
OUT_DIR="data/03_data_analysis/09_suppa2"
SHORT_DIR="data/02_build_shortread_transcriptome/05_edgeR"
LONG_DIR="data/01_build_longread_transcriptome/16_edgeR_filterexp"
mkdir -p "$OUT_DIR"

EVENTS=("A3" "A5" "AF" "AL" "MX" "RI" "SE")

# Keep only the columns of one sex (header has sample names only; rows are ID + 12 samples)
# Usage: split_sex <pattern> <in> <out>
split_sex() {
    awk -v pat="$1" 'BEGIN{FS="\t"; OFS="\t"} {
        if (NR==1) {
            found=0; for(i=1; i<=NF; i++) if($i ~ pat) { printf (found++ ? OFS : "") $i; h[i]=1 }
            printf "\n"
        } else {
            printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
        }
    }' "$2" > "$3"
    # SUPPA2 wants the header shifted by one tab
    perl -i -pe 's/^/\t/ if $. == 1' "$3"
}

# Run the whole analysis for one tissue
# Usage: run_tissue <tissue> <short_csv> <long_csv> <female_pattern> <male_pattern> <female_label> <male_label>
run_tissue() {
    local TISSUE=$1 SHORT_CSV=$2 LONG_CSV=$3 F_PAT=$4 M_PAT=$5 F_LAB=$6 M_LAB=$7
    local SUMMARY_FILE="${OUT_DIR}/events_summary_${TISSUE}.tsv"

    if [[ -f "$SUMMARY_FILE" ]]; then
        echo "$TISSUE already done ($SUMMARY_FILE exists), skipping"
        return
    fi

    echo "=========================================="
    echo "TISSUE: $TISSUE"
    echo "=========================================="

    # --- (1) TRANSFORM TPM MATRICES INTO SUPPA2 FORMAT ---

    # Short reads: first column is the row name; find the "transcript" column for the ID
    sed 's/"//g' "$SHORT_CSV" | \
    awk 'BEGIN{FS=","; OFS="\t"}
    {
        if (NR == 1) {
            for (i=1; i<=NF; i++) {
                if ($i == "transcript") col_idx = i
            }
            # Header: sample names only (columns 2-13)
            print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
        } else {
            # ID, then the 12 samples
            print $col_idx, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
        }
    }' > "${OUT_DIR}/short_tpm_${TISSUE}.tsv"

    # Long reads: samples are columns 1-12; ID is the second-to-last column
    sed 's/"//g' "$LONG_CSV" | \
    awk 'BEGIN{FS=","; OFS="\t"}
    {
        if (NR == 1) {
            # Header: sample names only
            printf "%s", $1
            for(i=2; i<=12; i++) printf "\t%s", $i
            printf "\n"
        } else {
            printf "%s", $(NF-1)
            for(i=1; i<=12; i++) printf "\t%s", $i
            printf "\n"
        }
    }' > "${OUT_DIR}/long_tpm_${TISSUE}.tsv"

    for LABEL in "short" "long"; do
        if [[ "$LABEL" == "short" ]]; then
            GTF="data/GTF/06_shortreadpipe_filtered_by_exp_short.gtf"
        else
            GTF="data/GTF/04_filtered_by_exp_long.gtf"
        fi
        TPM="${OUT_DIR}/${LABEL}_tpm_${TISSUE}.tsv"
        SUB_DIR="${OUT_DIR}/${LABEL}_${TISSUE}"
        mkdir -p "$SUB_DIR/events" "$SUB_DIR/psi" "$SUB_DIR/diff"

        echo "------------------------------------------"
        echo "$TISSUE: $LABEL reads"
        echo "------------------------------------------"

        # --- (2) GENERATE EVENTS AND CALCULATE PSI ---
        # Long-read GTFs often have more novel isoforms
        python3 $SUPPA_PATH generateEvents \
            -i "$GTF" \
            -o "${SUB_DIR}/events/${LABEL}_events" \
            -f ioe -e SE SS MX RI FL \
            -b V \
            -t 0

        for ioe_file in "${SUB_DIR}"/events/${LABEL}_events_*.ioe; do
            # Event type from the file name (e.g. short_events_SE_variable_0.ioe -> SE)
            EVENT_TYPE=$(basename "$ioe_file" | sed -E "s/${LABEL}_events_([A-Z0-9]+)_.*/\1/")
            echo "Calculating PSI for $LABEL - $EVENT_TYPE"
            python3 $SUPPA_PATH psiPerEvent \
                -i "$ioe_file" \
                -e "$TPM" \
                -o "${SUB_DIR}/psi/${LABEL}_${EVENT_TYPE}"
        done

        # --- (3) SPLIT BY SEX AND RUN DIFFSPLICE ---
        split_sex "$F_PAT" "$TPM" "${SUB_DIR}/${LABEL}_${F_LAB}.tpm"
        split_sex "$M_PAT" "$TPM" "${SUB_DIR}/${LABEL}_${M_LAB}.tpm"

        for EVENT in "${EVENTS[@]}"; do
            IOE_FILE="${SUB_DIR}/events/${LABEL}_events_${EVENT}_variable_0.ioe"
            PSI_INPUT="${SUB_DIR}/psi/${LABEL}_${EVENT}.psi"

            if [[ -f "$IOE_FILE" ]]; then
                echo "--- Working on $EVENT ($LABEL) ---"
                split_sex "$F_PAT" "$PSI_INPUT" "${SUB_DIR}/psi/${LABEL}_${F_LAB}_${EVENT}.psi"
                split_sex "$M_PAT" "$PSI_INPUT" "${SUB_DIR}/psi/${LABEL}_${M_LAB}_${EVENT}.psi"

                python3 $SUPPA_PATH diffSplice \
                    -m empirical \
                    -i "$IOE_FILE" \
                    -p "${SUB_DIR}/psi/${LABEL}_${F_LAB}_${EVENT}.psi" "${SUB_DIR}/psi/${LABEL}_${M_LAB}_${EVENT}.psi" \
                    -e "${SUB_DIR}/${LABEL}_${F_LAB}.tpm" "${SUB_DIR}/${LABEL}_${M_LAB}.tpm" \
                    --lower-bound 0 \
                    -nan 0.5 \
                    -gc \
                    -o "${SUB_DIR}/diff/${LABEL}_${EVENT}"
            else
                echo "Warning: Could not find $IOE_FILE"
            fi
        done

        # Rename .temp.0 files to proper .dpsi files (SUPPA output quirk)
        for f in "${SUB_DIR}"/diff/*.dpsi.temp.0; do
            if [[ -e "$f" ]]; then mv "$f" "${f%.temp.0}"; fi
        done
    done

    # --- (4) GATHER EVENTS FOR THIS TISSUE ---
    echo -e "mode\tEvent_Type\tEvent_ID\tdPSI\tp_value" > "$SUMMARY_FILE"
    for LABEL in "short" "long"; do
        for EVENT in "${EVENTS[@]}"; do
            FILE="${OUT_DIR}/${LABEL}_${TISSUE}/diff/${LABEL}_${EVENT}.dpsi"
            if [[ -f "$FILE" ]]; then
                # NR>1 skips the header line of the individual file
                awk -v m="$LABEL" -v e="$EVENT" 'NR>1 {print m"\t"e"\t"$0}' "$FILE" >> "$SUMMARY_FILE"
            else
                echo "Warning: $FILE not found, skipping..."
            fi
        done
    done

    echo "Done! Summary for $TISSUE: $SUMMARY_FILE ($(wc -l < "$SUMMARY_FILE") lines including header)"
}

run_tissue "gonad" \
    "${SHORT_DIR}/01_gonadtran_TMMTPMmatrix_short.csv" \
    "${LONG_DIR}/01_gonadtran_TMMTPMmatrix.csv" \
    "^RFG" "^RMG" "FGP" "MGP"

run_tissue "tail" \
    "${SHORT_DIR}/02_tailtran_TMMTPMmatrix_short.csv" \
    "${LONG_DIR}/02_tailtran_TMMTPMmatrix.csv" \
    "^RFT" "^RMT" "FTP" "MTP"

# ---------------------------------------------------------
# JOIN TAIL AND GONAD RESULTS
# ---------------------------------------------------------
echo "------------------------------------------"
echo "COMBINING GONAD AND TAIL RESULTS"
echo "------------------------------------------"

FINAL_COMBINED="${OUT_DIR}/events_summary_combined.tsv"
echo -e "tissue\tmode\tEvent_Type\tEvent_ID\tdPSI\tp_value" > "$FINAL_COMBINED"
for TISSUE in "gonad" "tail"; do
    awk -v t="$TISSUE" 'NR>1 {print t"\t"$0}' "${OUT_DIR}/events_summary_${TISSUE}.tsv" >> "$FINAL_COMBINED"
done

echo "Combined results saved to: $FINAL_COMBINED ($(wc -l < "$FINAL_COMBINED") lines including header)"
