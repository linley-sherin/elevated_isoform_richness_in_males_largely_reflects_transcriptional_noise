## Step 10: match long-read and short-read isoforms by intron chain (1-to-1 best fit)
## Last modified by Linley Sherin 2026-10-23
#
# Usage (from the project folder):
#   python3 03_data_analysis/10_match_transcripts.py 0
# where '0' is the allowance (bp) for fuzzy junctions (0 = exact; used for the paper)
#
# Reads:  data/GTF/06_shortreadpipe_filtered_by_exp_short.gtf   (short-read step 05)
#         data/GTF/04_filtered_by_exp_long.gtf                  (long-read step 16)
# Writes: data/03_data_analysis/10_match_transcripts/
#           <tolerance>_matched_ids.csv    gene_id, long_read_id, short_read_id, status; input to step 11
#           <tolerance>_venn_results.csv   shared / long-read unique / short-read unique counts

import sys
import os
import csv

def get_introns(gtf_file):
    """Parses GTF and returns a dict of {(gene_id, transcript_id): intron_chain_tuple}"""
    transcripts = {}
    print(f"Reading {gtf_file}...")
    with open(gtf_file, 'r') as f:
        for line in f:
            if line[0] == '#' or 'exon' not in line:
                continue
            parts = line.split('\t')
            if parts[2] != 'exon':
                continue
            
            # Attribute parsing
            attr = parts[8]
            try:
                gid = attr.split('gene_id "')[1].split('"')[0]
                tid = attr.split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue
            
            key = (gid, tid)
            if key not in transcripts:
                transcripts[key] = []
            transcripts[key].append((int(parts[3]), int(parts[4])))
    
    intron_chains = {}
    for (gid, tid), exons in transcripts.items():
        if len(exons) > 1:
            exons.sort()
            # Tuple conversion for O(1) comparison and hashing
            chain = tuple((exons[i][1], exons[i+1][0]) for i in range(len(exons) - 1))
            intron_chains[(gid, tid)] = chain
    return intron_chains

def is_fuzzy_match(chain1, chain2, tol):
    """Returns True if all junctions are within the specified tolerance."""
    for i in range(len(chain1)):
        if abs(chain1[i][0] - chain2[i][0]) > tol or abs(chain1[i][1] - chain2[i][1]) > tol:
            return False
    return True

def get_match_score(chain1, chain2):
    """Calculates total offset between two chains. Lower score = better fit."""
    total_offset = 0
    for i in range(len(chain1)):
        total_offset += abs(chain1[i][0] - chain2[i][0])
        total_offset += abs(chain1[i][1] - chain2[i][1])
    return total_offset

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 03_data_analysis/10_match_transcripts.py <tolerance>")
        print("Example: python3 03_data_analysis/10_match_transcripts.py 0")
        sys.exit(1)

    tolerance = int(sys.argv[1])

    # Paths are relative to the project root (one folder up from this script)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sr_gtf = os.path.join(root, "data", "GTF", "06_shortreadpipe_filtered_by_exp_short.gtf")
    lr_gtf = os.path.join(root, "data", "GTF", "04_filtered_by_exp_long.gtf")
    out_dir = os.path.join(root, "data", "03_data_analysis", "10_match_transcripts")
    os.makedirs(out_dir, exist_ok=True)

    # 1. Load Data
    sr_data = get_introns(sr_gtf)
    lr_data = get_introns(lr_gtf)

    # Index SR data by gene for faster lookup
    sr_by_gene = {}
    for (gid, tid), chain in sr_data.items():
        sr_by_gene.setdefault(gid, []).append((tid, chain))

    # 2. Matching Logic (Best-Fit / 1-to-1)
    match_map = {}      # {lr_tid: sr_tid}
    shared_sr_keys = set() # {(sr_gid, sr_tid)}
    used_sr_tids = set()   # {sr_tid} to prevent multiple matches

    print(f"Matching within genes (Tolerance: {tolerance}bp, Best-Fit)...")
    count = 0
    
    # Iterate through Long Reads
    for (lr_gid, lr_tid), lr_chain in lr_data.items():
        count += 1
        if count % 10000 == 0:
            print(f"  Processed {count} LR transcripts...")

        if lr_gid in sr_by_gene:
            best_sr_tid = None
            min_offset = float('inf')

            # Look for the best match among Short Reads in the same gene
            for sr_tid, sr_chain in sr_by_gene[lr_gid]:
                # Skip if this SR transcript was already "claimed" by another LR transcript
                if sr_tid in used_sr_tids:
                    continue

                # Must have same number of introns
                if len(lr_chain) != len(sr_chain):
                    continue
                
                # Check if within fuzzy tolerance
                if is_fuzzy_match(lr_chain, sr_chain, tolerance):
                    offset = get_match_score(lr_chain, sr_chain)
                    
                    # If this is the closest match so far, track it
                    if offset < min_offset:
                        min_offset = offset
                        best_sr_tid = sr_tid
            
            # If a best match was found, finalize it
            if best_sr_tid:
                match_map[lr_tid] = best_sr_tid
                shared_sr_keys.add((lr_gid, best_sr_tid))
                used_sr_tids.add(best_sr_tid)

    # 3. Calculate Totals
    total_lr = len(lr_data)
    total_sr = len(sr_data)
    intersect = len(match_map)
    lr_unique = total_lr - intersect
    sr_unique = total_sr - len(shared_sr_keys)
    
    mapping_file = os.path.join(out_dir, f"{tolerance}_matched_ids.csv")
    venn_file = os.path.join(out_dir, f"{tolerance}_venn_results.csv")

    # 4. Output Comprehensive Mapping (including gene_id)
    with open(mapping_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['gene_id', 'long_read_id', 'short_read_id', 'status'])
        
        # Write Matched and LR-Unique
        for (lr_gid, lr_tid) in lr_data.keys():
            if lr_tid in match_map:
                writer.writerow([lr_gid, lr_tid, match_map[lr_tid], 'matched'])
            else:
                writer.writerow([lr_gid, lr_tid, 'NA', 'lr_unique'])
        
        # Write SR-Unique
        for (sr_gid, sr_tid) in sr_data.keys():
            if (sr_gid, sr_tid) not in shared_sr_keys:
                writer.writerow([sr_gid, 'NA', sr_tid, 'sr_unique'])

    # 5. Output Venn Summary
    with open(venn_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Count'])
        writer.writerow(['Intersection', intersect])
        writer.writerow(['Long-Read_Unique', lr_unique])
        writer.writerow(['Short-Read_Unique', sr_unique])

    print("\n" + "="*55)
    print(f"RESULTS (Tolerance: {tolerance}bp | 1-to-1 Best Fit)")
    print(f"Shared: {intersect} | LR Unique: {lr_unique} | SR Unique: {sr_unique}")
    print("="*55)
    print(f"Mapping saved to: {mapping_file}")
    print(f"Venn summary saved to: {venn_file}")
