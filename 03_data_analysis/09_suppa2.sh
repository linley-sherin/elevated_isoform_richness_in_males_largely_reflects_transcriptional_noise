#!/bin/bash


# ---------------------------------------------------------
# ---------------------------------------------------------
# GONAD
# ---------------------------------------------------------
# ---------------------------------------------------------
# --- TRANFORM TPM MATRICES INTO PROPER FORMAT ---

# SHORT GONAD
sed 's/"//g' data/02_build_shortread_transcriptome/05_edgeR/01_gonadtran_TMMTPMmatrix_short.csv | \
awk 'BEGIN{FS=","; OFS="\t"} 
{
    if (NR == 1) {
        # 1. Find which column index is "transcript"
        for (i=1; i<=NF; i++) {
            if ($i == "transcript") col_idx = i
        }
        # 2. Print ONLY the sample names for the header (matching your ENST example)
        print "1FGP", "1MGP", "2FGP", "2MGP", "3FGP", "3MGP", "4FGP", "4MGP", "5FGP", "5MGP", "6FGP", "6MGP"
    } else {
        # 3. Print the ID (from the index we found) followed by the sample data (cols 2-13)
        # This ensures we have 1 (ID) + 12 (Samples) = 13 total columns
        print $col_idx, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
    }
}' > ./15_suppa2/short_tpm_gonad.tsv

# LONG GONAD
sed 's/"//g' data/01_build_longread_transcriptome/16_edgeR_filterexp/01_gonadtran_TMMTPMmatrix.csv | \
awk 'BEGIN{FS=","; OFS="\t"} 
{
    if (NR == 1) {
        # HEADER: 12 sample names
        print "1FGP", "1MGP", "2FGP", "2MGP", "3FGP", "3MGP", "4FGP", "4MGP", "5FGP", "5MGP", "6FGP", "6MGP"
    } else {
        # DATA: 
        # 1. Grab ID from the second-to-last column $(NF-1)
        # 2. Print that ID, then numeric samples 1 through 12
        
        printf "%s", $(NF-1);
        for(i=1; i<=12; i++) {
            printf "\t%s", $i;
        }
        printf "\n";
    }
}' > ./15_suppa2/long_tpm_gonad.tsv

# --- GLOBAL CONFIG ---
SUPPA_PATH="tools/SUPPA/suppa.py"
OUT_DIR="./15_suppa2/"
mkdir -p $OUT_DIR

# --- DATASET DEFINITIONS ---
# Format: "Label|GTF_Path|TPM_Path"
DATASETS=(
    "short|data/GTF/06_shortreadpipe_filtered_by_exp_short.gtf|./15_suppa2/short_tpm_gonad.tsv"
    "long|data/GTF/04_filtered_by_exp_long.gtf|./15_suppa2/long_tpm_gonad.tsv"
)

# --- PROCESSING LOOP ---
for data in "${DATASETS[@]}"; do
    IFS="|" read -r LABEL GTF TPM <<< "$data"
    
    echo "------------------------------------------"
    echo "STARTING ANALYSIS FOR: $LABEL"
    echo "------------------------------------------"

    # Create specific subdirectories
    SUB_DIR="${OUT_DIR}/${LABEL}"
    mkdir -p "$SUB_DIR/events" "$SUB_DIR/psi"

    # 1. Generate Events for this specific GTF
    # Note: Long-read GTFs often have more novel isoforms
    python3 $SUPPA_PATH generateEvents \
        -i "$GTF" \
        -o "${SUB_DIR}/events/${LABEL}_events" \
        -f ioe -e SE SS MX RI FL \
        -b V \
        -t 0

	# 2. Calculate PSI Values
    for ioe_file in ${SUB_DIR}/events/${LABEL}_events_*.ioe; do
        # This regex strips everything from "_variable" onwards to get just "SE", "A3", etc.
        EVENT_TYPE=$(basename "$ioe_file" | sed -E "s/${LABEL}_events_([A-Z0-9]+)_.*/\1/")
        
        echo "Calculating PSI for $LABEL - $EVENT_TYPE"
        
        # FIX: We explicitly name the output so it matches the downstream expectation: short_SE.psi
        python3 $SUPPA_PATH psiPerEvent \
            -i "$ioe_file" \
            -e "$TPM" \
            -o "${SUB_DIR}/psi/${LABEL}_${EVENT_TYPE}"
    done

    echo "Finished processing $LABEL reads."
done

echo "=========================================="
echo "All datasets processed! Results are in: $OUT_DIR"

# ---------------------------------------------------------
# DYNAMIC SPLITTING AND DIFFSPLICE
# ---------------------------------------------------------

for MODE in "short" "long"; do
    echo "=========================================="
    echo "PROCESSING $MODE DATA"
    echo "=========================================="
    
    mkdir -p ./15_suppa2/${MODE}/diff
    TPM_FILE="./15_suppa2/${MODE}_tpm_gonad.tsv"

    # 1. SPLIT TPMs
    awk 'BEGIN{FS="\t"; OFS="\t"} {
        if (NR==1) {
            found=0; for(i=1; i<=NF; i++) if($i ~ /FGP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
            printf "\n"
        } else {
            printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
        }
    }' $TPM_FILE > ./15_suppa2/${MODE}/${MODE}_FGP.tpm

    awk 'BEGIN{FS="\t"; OFS="\t"} {
        if (NR==1) {
            found=0; for(i=1; i<=NF; i++) if($i ~ /MGP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
            printf "\n"
        } else {
            printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
        }
    }' $TPM_FILE > ./15_suppa2/${MODE}/${MODE}_MGP.tpm

    perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/${MODE}_FGP.tpm
    perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/${MODE}_MGP.tpm

    for EVENT in "A3" "A5" "AF" "AL" "MX" "RI" "SE"; do
        # UPDATED PATH: Matches your 'ls' output (_variable_0.ioe)
        IOE_FILE="./15_suppa2/${MODE}/events/${MODE}_events_${EVENT}_variable_0.ioe"
        PSI_INPUT="./15_suppa2/${MODE}/psi/${MODE}_${EVENT}.psi"

        if [[ -f "$IOE_FILE" ]]; then
            echo "--- Working on $EVENT ($MODE) ---"
            
            # 2. SPLIT PSIs
            awk 'BEGIN{FS="\t"; OFS="\t"} {
                if (NR==1) {
                    found=0; for(i=1; i<=NF; i++) if($i ~ /FGP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
                    printf "\n"
                } else {
                    printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
                }
            }' $PSI_INPUT > ./15_suppa2/${MODE}/psi/${MODE}_FGP_${EVENT}.psi

            awk 'BEGIN{FS="\t"; OFS="\t"} {
                if (NR==1) {
                    found=0; for(i=1; i<=NF; i++) if($i ~ /MGP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
                    printf "\n"
                } else {
                    printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
                }
            }' $PSI_INPUT > ./15_suppa2/${MODE}/psi/${MODE}_MGP_${EVENT}.psi

            perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/psi/${MODE}_FGP_${EVENT}.psi
            perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/psi/${MODE}_MGP_${EVENT}.psi

            # 3. RUN DIFFSPLICE
            python3 tools/SUPPA/suppa.py diffSplice \
                -m empirical \
                -i "$IOE_FILE" \
                -p ./15_suppa2/${MODE}/psi/${MODE}_FGP_${EVENT}.psi ./15_suppa2/${MODE}/psi/${MODE}_MGP_${EVENT}.psi \
                -e ./15_suppa2/${MODE}/${MODE}_FGP.tpm ./15_suppa2/${MODE}/${MODE}_MGP.tpm \
                --lower-bound 0 \
                -nan 0.5 \
                -gc \
                -o ./15_suppa2/${MODE}/diff/${MODE}_${EVENT}
        else
            echo "Warning: Could not find $IOE_FILE"
        fi
    done
done
# ---------------------------------------------------------
# SUMMARY: Gathering Events
# ---------------------------------------------------------
echo "------------------------------------------"
echo "GENERATING RESULTS SUMMARY"
echo "------------------------------------------"

# Rename all .temp.0 files to proper .dpsi files (SUPPA output quirk)
for f in ./15_suppa2/*/diff/*.dpsi.temp.0; do
    [ -e "$f" ] && mv "$f" "${f%.temp.0}"
done

SUMMARY_FILE="./15_suppa2/events_summary.tsv"
echo -e "mode\tEvent_Type\tEvent_ID\tdPSI\tp_value" > $SUMMARY_FILE

for MODE in "short" "long"; do
    for EVENT in "A3" "A5" "AF" "AL" "MX" "RI" "SE"; do
        FILE="./15_suppa2/${MODE}/diff/${MODE}_${EVENT}.dpsi"
        
        if [[ -f "$FILE" ]]; then
            # NR>1 skips the header line of the individual file
            # print m, e, and the whole line ($0)
            awk -v m="$MODE" -v e="$EVENT" 'NR>1 {print m"\t"e"\t"$0}' "$FILE" >> $SUMMARY_FILE
        else
            echo "Warning: $FILE not found, skipping..."
        fi
    done
done

echo "Done! Summary saved in: $SUMMARY_FILE"
echo "Total lines (including header):"
wc -l $SUMMARY_FILE

# ---------------------------------------------------------
# FINAL RENAMING (The "_gonad" part)
# ---------------------------------------------------------
echo "Renaming directories and summary file..."

# Rename the final summary
mv ./15_suppa2/events_summary.tsv ./15_suppa2/events_summary_gonad.tsv

# Rename the main results folders
[ -d "./15_suppa2/short" ] && mv ./15_suppa2/short ./15_suppa2/short_gonad
[ -d "./15_suppa2/long" ] && mv ./15_suppa2/long ./15_suppa2/long_gonad

echo "Done! Final summary: ./15_suppa2/events_summary_gonad.tsv"


# ---------------------------------------------------------
# ---------------------------------------------------------
# TAIL
# ---------------------------------------------------------
# ---------------------------------------------------------

# --- TRANFORM TPM MATRICES INTO PROPER FORMAT ---

# SHORT tail
sed 's/"//g' data/02_build_shortread_transcriptome/05_edgeR/02_tailtran_TMMTPMmatrix_short.csv | \
awk 'BEGIN{FS=","; OFS="\t"} 
{
    if (NR == 1) {
        # 1. Find which column index is "transcript"
        for (i=1; i<=NF; i++) {
            if ($i == "transcript") col_idx = i
        }
        # 2. Print ONLY the sample names for the header (matching your ENST example)
        print "1FTP", "1MTP", "2FTP", "2MTP", "3FTP", "3MTP", "4FTP", "4MTP", "5FTP", "5MTP", "6FTP", "6MTP"
    } else {
        # 3. Print the ID (from the index we found) followed by the sample data (cols 2-13)
        # This ensures we have 1 (ID) + 12 (Samples) = 13 total columns
        print $col_idx, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13
    }
}' > ./15_suppa2/short_tpm_tail.tsv

# LONG TAIL
sed 's/"//g' data/01_build_longread_transcriptome/16_edgeR_filterexp/02_tailtran_TMMTPMmatrix.csv | \
awk 'BEGIN{FS=","; OFS="\t"} 
{
    if (NR == 1) {
        # HEADER: Exactly 12 sample names.
        print "1FTP", "1MTP", "2FTP", "2MTP", "3FTP", "3MTP", "4FTP", "4MTP", "5FTP", "5MTP", "6FTP", "6MTP"
    } else {
        # DATA: 13 columns
        printf "%s", $1;
        for(i=2; i<=13; i++) {
            printf "\t%s", $i;
        }
        printf "\n";
    }
}' > ./15_suppa2/long_tpm_tail.tsv

# --- GLOBAL CONFIG ---
SUPPA_PATH="tools/SUPPA/suppa.py"
OUT_DIR="./15_suppa2/"
mkdir -p $OUT_DIR

# --- DATASET DEFINITIONS ---
# Format: "Label|GTF_Path|TPM_Path"
DATASETS=(
    "short|data/GTF/06_shortreadpipe_filtered_by_exp_short.gtf|./15_suppa2/short_tpm_tail.tsv"
    "long|data/GTF/04_filtered_by_exp_long.gtf|./15_suppa2/long_tpm_tail.tsv"
)

# --- PROCESSING LOOP ---
for data in "${DATASETS[@]}"; do
    IFS="|" read -r LABEL GTF TPM <<< "$data"
    
    echo "------------------------------------------"
    echo "STARTING ANALYSIS FOR: $LABEL"
    echo "------------------------------------------"

    # Create specific subdirectories
    SUB_DIR="${OUT_DIR}/${LABEL}"
    mkdir -p "$SUB_DIR/events" "$SUB_DIR/psi"

    # 1. Generate Events for this specific GTF
    # Note: Long-read GTFs often have more novel isoforms
    python3 $SUPPA_PATH generateEvents \
        -i "$GTF" \
        -o "${SUB_DIR}/events/${LABEL}_events" \
        -f ioe -e SE SS MX RI FL \
        -b V \
        -t 0

	# 2. Calculate PSI Values
    for ioe_file in ${SUB_DIR}/events/${LABEL}_events_*.ioe; do
        # This regex strips everything from "_variable" onwards to get just "SE", "A3", etc.
        EVENT_TYPE=$(basename "$ioe_file" | sed -E "s/${LABEL}_events_([A-Z0-9]+)_.*/\1/")
        
        echo "Calculating PSI for $LABEL - $EVENT_TYPE"
        
        # FIX: We explicitly name the output so it matches the downstream expectation: short_SE.psi
        python3 $SUPPA_PATH psiPerEvent \
            -i "$ioe_file" \
            -e "$TPM" \
            -o "${SUB_DIR}/psi/${LABEL}_${EVENT_TYPE}"
    done

    echo "Finished processing $LABEL reads."
done

echo "=========================================="
echo "All datasets processed! Results are in: $OUT_DIR"

# ---------------------------------------------------------
# DYNAMIC SPLITTING AND DIFFSPLICE
# ---------------------------------------------------------

for MODE in "short" "long"; do
    echo "=========================================="
    echo "PROCESSING $MODE DATA"
    echo "=========================================="
    
    mkdir -p ./15_suppa2/${MODE}/diff
    TPM_FILE="./15_suppa2/${MODE}_tpm_tail.tsv"

    # 1. SPLIT TPMs
    awk 'BEGIN{FS="\t"; OFS="\t"} {
        if (NR==1) {
            found=0; for(i=1; i<=NF; i++) if($i ~ /FTP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
            printf "\n"
        } else {
            printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
        }
    }' $TPM_FILE > ./15_suppa2/${MODE}/${MODE}_FTP.tpm

    awk 'BEGIN{FS="\t"; OFS="\t"} {
        if (NR==1) {
            found=0; for(i=1; i<=NF; i++) if($i ~ /MTP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
            printf "\n"
        } else {
            printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
        }
    }' $TPM_FILE > ./15_suppa2/${MODE}/${MODE}_MTP.tpm

    perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/${MODE}_FTP.tpm
    perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/${MODE}_MTP.tpm

    for EVENT in "A3" "A5" "AF" "AL" "MX" "RI" "SE"; do
        # UPDATED PATH: Matches your 'ls' output (_variable_0.ioe)
        IOE_FILE="./15_suppa2/${MODE}/events/${MODE}_events_${EVENT}_variable_0.ioe"
        PSI_INPUT="./15_suppa2/${MODE}/psi/${MODE}_${EVENT}.psi"

        if [[ -f "$IOE_FILE" ]]; then
            echo "--- Working on $EVENT ($MODE) ---"
            
            # 2. SPLIT PSIs
            awk 'BEGIN{FS="\t"; OFS="\t"} {
                if (NR==1) {
                    found=0; for(i=1; i<=NF; i++) if($i ~ /FTP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
                    printf "\n"
                } else {
                    printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
                }
            }' $PSI_INPUT > ./15_suppa2/${MODE}/psi/${MODE}_FTP_${EVENT}.psi

            awk 'BEGIN{FS="\t"; OFS="\t"} {
                if (NR==1) {
                    found=0; for(i=1; i<=NF; i++) if($i ~ /MTP/) { printf (found++ ? OFS : "") $i; h[i]=1 }
                    printf "\n"
                } else {
                    printf $1; for(i=1; i<=12; i++) if(h[i]) printf "\t"$(i+1); printf "\n"
                }
            }' $PSI_INPUT > ./15_suppa2/${MODE}/psi/${MODE}_MTP_${EVENT}.psi

            perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/psi/${MODE}_FTP_${EVENT}.psi
            perl -i -pe 's/^/\t/ if $. == 1' ./15_suppa2/${MODE}/psi/${MODE}_MTP_${EVENT}.psi

            # 3. RUN DIFFSPLICE
            python3 tools/SUPPA/suppa.py diffSplice \
                -m empirical \
                -i "$IOE_FILE" \
                -p ./15_suppa2/${MODE}/psi/${MODE}_FTP_${EVENT}.psi ./15_suppa2/${MODE}/psi/${MODE}_MTP_${EVENT}.psi \
                -e ./15_suppa2/${MODE}/${MODE}_FTP.tpm ./15_suppa2/${MODE}/${MODE}_MTP.tpm \
                --lower-bound 0 \
                -nan 0.5 \
                -gc \
                -o ./15_suppa2/${MODE}/diff/${MODE}_${EVENT}
        else
            echo "Warning: Could not find $IOE_FILE"
        fi
    done
done

# ---------------------------------------------------------
# SUMMARY: Gathering Events
# ---------------------------------------------------------
echo "------------------------------------------"
echo "GENERATING RESULTS SUMMARY"
echo "------------------------------------------"

# Rename all .temp.0 files to proper .dpsi files (SUPPA output quirk)
for f in ./15_suppa2/*/diff/*.dpsi.temp.0; do
    [ -e "$f" ] && mv "$f" "${f%.temp.0}"
done

SUMMARY_FILE="./15_suppa2/events_summary.tsv"
echo -e "mode\tEvent_Type\tEvent_ID\tdPSI\tp_value" > $SUMMARY_FILE

for MODE in "short" "long"; do
    for EVENT in "A3" "A5" "AF" "AL" "MX" "RI" "SE"; do
        FILE="./15_suppa2/${MODE}/diff/${MODE}_${EVENT}.dpsi"
        
        if [[ -f "$FILE" ]]; then
            # NR>1 skips the header line of the individual file
            # print m, e, and the whole line ($0)
            awk -v m="$MODE" -v e="$EVENT" 'NR>1 {print m"\t"e"\t"$0}' "$FILE" >> $SUMMARY_FILE
        else
            echo "Warning: $FILE not found, skipping..."
        fi
    done
done

echo "Done! Summary saved in: $SUMMARY_FILE"
echo "Total lines (including header):"
wc -l $SUMMARY_FILE

# ---------------------------------------------------------
# FINAL RENAMING (The "_tail" part)
# ---------------------------------------------------------
echo "Renaming directories and summary file..."

# Rename the final summary
mv ./15_suppa2/events_summary.tsv ./15_suppa2/events_summary_tail.tsv

# Rename the main results folders
[ -d "./15_suppa2/short" ] && mv ./15_suppa2/short ./15_suppa2/short_tail
[ -d "./15_suppa2/long" ] && mv ./15_suppa2/long ./15_suppa2/long_tail

echo "Done! Final summary: ./15_suppa2/events_summary_tail.tsv"


# ---------------------------------------------------------
# ---------------------------------------------------------
# JOIN TAIL AND GONAD RESULTS
# ---------------------------------------------------------
# ---------------------------------------------------------


echo "------------------------------------------"
echo "COMBINING GONAD AND TAIL RESULTS"
echo "------------------------------------------"

FINAL_COMBINED="./15_suppa2/events_summary_combined.tsv"

# 1. Create the header for the new combined file
echo -e "tissue\tmode\tEvent_Type\tEvent_ID\tdPSI\tp_value" > $FINAL_COMBINED

# 2. Add Gonad data (NR>1 skips the header of the source file)
if [[ -f "./15_suppa2/events_summary_gonad.tsv" ]]; then
    awk 'NR>1 {print "gonad\t"$0}' ./15_suppa2/events_summary_gonad.tsv >> $FINAL_COMBINED
else
    echo "Error: gonad summary file not found!"
fi

# 3. Add Tail data
if [[ -f "./15_suppa2/events_summary_tail.tsv" ]]; then
    awk 'NR>1 {print "tail\t"$0}' ./15_suppa2/events_summary_tail.tsv >> $FINAL_COMBINED
else
    echo "Error: tail summary file not found!"
fi

echo "Success! Combined results saved to: $FINAL_COMBINED"
echo "Total events found:"
wc -l $FINAL_COMBINED