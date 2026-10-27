## Step 14: remove redundant and degraded isoforms from the SQANTI-filtered transcriptome
## Last modified by Linley Sherin 2026-10-27
#
# Despite attempting to collapse redundancy and 5' end, the gtf output from
# SQANTI3 still has both. So here I am further filtering my gtf/transcriptome
# before quantifying to remove redundancy and potentially degraded transcripts
#
# Runs locally.
# (1) Load
# (2) Create SJ whitelist so that LRT and SRT are comparable
# (3) Correct exon-structure redundancy
# (4) Snap 5' and 3' ends to consensus boundaries
# (5) Correct 5' collapse (degradation)
# (x) Collapse ORF redundancy - DON'T DO
# (x) Re-cluster within pb_ids to determine if there are multiple genes
# (6) Switch gene names to reference names using SQANTI output
# (7) Remove non-coding genes
# (8) Remove non-autosomal genes
# (9) Save results
# (10) Filter fasta to match
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 01_build_longread_transcriptome/14_remove_redundancy.R
#
# Reads:  data/GTF/02_SQANTI_filtered.gtf
#         data/01_build_longread_transcriptome/06_mapcollapse/collapsed.abundance.txt
#         data/01_build_longread_transcriptome/09_star/*SJ.out.tab
#         data/01_build_longread_transcriptome/11_sqanti_qc/output/SQANTI_output_junctions.txt
#         data/01_build_longread_transcriptome/12_sqanti_filter/
#           SQANTI_filtered_RulesFilter_result_classification.txt, SQANTI_filtered.filtered.fasta
#         data/01_build_longread_transcriptome/13_transdecoder/parsed_results.csv
# Writes: data/01_build_longread_transcriptome/14_remove_redundancy/
#           junction_whitelist.csv
#           nonredundant_transcripts.fasta   input to step 15
#           rosetta_stone.csv                old -> new gene/transcript IDs, input to step 16
#           03_nonredundant_transcriptome.gtf
#         data/GTF/03_nonredundant_transcriptome.gtf
#
# Overwrites its outputs every time it runs.

library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(rtracklayer)
library(stringr)
library(Biostrings)
library(purrr)

################################################################################
### (1) Load
################################################################################

# Read filtered gtf
gtf_sqanti <- rtracklayer::import(here::here("data",
                                             "GTF",
                                             "02_SQANTI_filtered.gtf"))

gtf_sqanti <- as.data.frame(gtf_sqanti) 

# Starting transcriptome info
cat("Before reducing redundancy, transcriptome contains", length(unique(gtf_sqanti$gene_id)), "unique genes.\n")
cat("Before reducing redundancy, transcriptome contains", length(unique(gtf_sqanti$transcript_id)), "unique transcripts.\n")

# Filter for transcript entries from gtf
transcripts <- gtf_sqanti[gtf_sqanti$type == "transcript", ]

# Read junction object using data.table (much faster for large SQANTI files)
junc <- data.table::fread(here::here("data",
                                     "01_build_longread_transcriptome",
                                     "11_sqanti_qc",
                                     "output",
                                     "SQANTI_output_junctions.txt"),
                          header = TRUE, 
                          sep = "\t", 
                          data.table = FALSE) # Returns a data frame to keep your script compatible

# Filter junction dataframe to keep only isoforms present in transcripts
junc_sqanti <- junc[junc$isoform %in% transcripts$transcript_id, ]

# Read expression data (FL transcripts from PacBio workflow)
fl_df <- data.table::fread(here::here("data",
                                "01_build_longread_transcriptome",
                                "06_mapcollapse",
                                "collapsed.abundance.txt"),
                     header = T,
                     skip = 4) %>%
  dplyr::select(transcript_id = pbid, fl_assoc) # keep only what you need from FL counts

# Read in orf data
orf_df <- read.csv(here::here("data",
                              "01_build_longread_transcriptome",
                              "13_transdecoder",
                              "parsed_results.csv"))

################################################################################
### (2) Create SJ whitelist so that LRT and SRT are comparable
################################################################################

# I found that even after SQANTI3 filtering on read coverage there
# are junctions not found in the STAR sj.out.tab files showing up
# in the LRT

# This is not really an issue except that for the SRT I need to remove all 
# non-STAR junctions so I am going to do the same for the LRT to keep them
# comparable

# SJ.out.tab directory
sj_dir <- here::here("data", "01_build_longread_transcriptome", "09_star")

# Get all SJ.out.tab files (non-recursive)
files <- list.files(
  path = sj_dir,
  pattern = "SJ\\.out\\.tab$",
  full.names = TRUE,
  recursive = FALSE
)

# Read all files and add sample column
sj <- files %>% 
  set_names() %>% 
  map_dfr(~ fread(.x) %>% 
            mutate(sample = basename(.x)), .id = "file")

# Name columns
colnames(sj)[2:10] <- c(
  "chr",
  "intron_start",
  "intron_end",
  "strand",
  "motif",
  "annotated",
  "unique_reads",
  "multi_reads",
  "max_overhang"
)

# Tidyverse pipeline to sum reads by junction and strand
sj_tidy <- sj %>%
  mutate(
    support = unique_reads + multi_reads,
    canonical = motif %in% c(1,2,3,5),
    strand_label = case_when(
      strand == 1 ~ "+",
      strand == 2 ~ "-",
      TRUE ~ "*"
    )
  ) %>%
  # Sum reads across all samples per junction
  group_by(chr, intron_start, intron_end, strand, motif, canonical, strand_label) %>%
  summarise(
    total_support = sum(support),
    .groups = "drop"
  ) %>%
  # Flag supported junctions (sum of reads >=10 and canonical)
  mutate(
    supported = total_support >= 10 & canonical
  )

# Subset to only supported junctions
sj_supported <- sj_tidy %>% filter(supported)

# Save
write_csv(sj_supported,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "14_remove_redundancy",
                     "junction_whitelist.csv"))

# Create a vector of supported junctions as chr-start-end strings
supported_juncs <- sj_supported %>%
  mutate(junction_coord = paste(chr, intron_start, intron_end, sep = "-")) %>%
  pull(junction_coord)

# Make supported junctions a tibble
supported_juncs_tbl <- sj_supported %>%
  transmute(junction_coord = paste(chr, intron_start, intron_end, strand_label, sep = "-"))

# Add junction_coord to each junction in junc_sqanti
junc_sqanti2 <- junc_sqanti %>%
  mutate(
    junction_coord = paste(chrom, genomic_start_coord, genomic_end_coord, strand, sep = "-")
  )

# Mark junctions that are supported
junc_sqanti2 <- junc_sqanti2 %>%
  left_join(supported_juncs_tbl %>% mutate(supported = TRUE),
            by = "junction_coord") %>%
  mutate(supported = replace_na(supported, FALSE))

# Identify removed isoforms
removed_isoforms <- junc_sqanti2 %>%
  group_by(isoform) %>%
  summarise(all_supported = all(supported), .groups = "drop") %>%
  filter(!all_supported) %>%
  pull(isoform)

# Keep only isoforms where all junctions are supported
iso_filtered <- junc_sqanti2 %>%
  group_by(isoform) %>%
  filter(all(supported)) %>%
  ungroup()

# Write removed isoforms to a log file
log_removed <- data.frame(
  transcript_id = removed_isoforms,
  reason = "non_supported_junctions",
  stringsAsFactors = FALSE
)

# Filter GTF to keep only isoforms with all junctions supported
gtf_step1 <- gtf_sqanti[gtf_sqanti$transcript_id %in% iso_filtered$isoform, ]

rm(removed_isoforms)
rm(supported_juncs)
rm(supported_juncs_tbl)
rm(junc_sqanti2)
rm(sj)
rm(sj_tidy)
rm(sj_supported)

################################################################################
### (3) Correct exon-structure redundancy
################################################################################

# Create a string of sorted junctions per isoform
junc_signature_df <- iso_filtered %>%
  mutate(junction_coord = paste0(genomic_start_coord, "-", genomic_end_coord)) %>%
  group_by(isoform, 
           chrom, 
           strand) %>%
  summarise(
    junction_signature = paste(sort(junction_coord), collapse = ";"),
    .groups = "drop"
  )

# Group by identical junction_signature
isoform_groups <- junc_signature_df %>%
  group_by(junction_signature) %>%
  summarise(
    isoforms = list(isoform),
    count = n(),
    .groups = "drop"
  )

### Of the collapsible isoforms, keep the most highly expressed isoform only

# Expand isoform groups into long format
isoform_groups <- isoform_groups %>%
  unnest(isoforms) %>%
  dplyr::rename(transcript_id = isoforms)

# Join expression data with isoform groups
isoform_groups_w_exp <- isoform_groups %>%
  left_join(fl_df, by = "transcript_id")

# Extract numeric part of isoform name for tie-breaking
isoform_groups_w_exp <- isoform_groups_w_exp %>%
  mutate(
    last_num = as.integer(str_extract(transcript_id, "(?<=\\.)\\d+$"))
  )

# Keep the most highly expressed isoform(s) per group
isoforms_to_keep <- isoform_groups_w_exp %>%
  group_by(junction_signature) %>%
  mutate(
    max_expr = max(fl_assoc, na.rm = TRUE) # Find max expression in group
  ) %>%
  filter(fl_assoc == max_expr) %>% # Only keep isoform that matches the max 
  slice_min(last_num, with_ties = FALSE) %>% # Selects the isoform with lowest # in id if there is a tie
  ungroup()

# Keep log of isoforms lost at this step and append to existing log_removed
log_removed <- bind_rows(
  log_removed,
  data.frame(
    transcript_id = setdiff(iso_filtered$isoform, isoforms_to_keep$transcript_id),
    reason = "exon_structure_redundancy",
    stringsAsFactors = FALSE
  )
)

# Filter gtf_step1 to keep only isoforms to keep
gtf_step2 <- gtf_step1[gtf_step1$transcript_id %in% isoforms_to_keep$transcript_id, ]

# Filter junction dataframe to keep only isoforms to keep
junc_step2 <- iso_filtered[iso_filtered$isoform %in% isoforms_to_keep$transcript_id, ]

# Remove unneeded objects to clean up R environment
rm(isoform_groups)
rm(isoform_groups_w_exp)
rm(junc_signature_df)
rm(isoforms_to_keep)

################################################################################
### (4) Snap 5' and 3' ends to consensus boundaries
################################################################################

# Because we can't trust the 5' and 3' ends completely and they cause
# differences downstream when doing event-level analysis we will 'snap' them all
# to match, using the most represented coordinate as the ground truth (the mode)

# Identify Terminal Exons
terminal_exons <- gtf_step2 %>%
  filter(type == "exon") %>%
  group_by(transcript_id) %>%
  mutate(
    exon_rank = if_else(strand == "+", rank(start), rank(-end)),
    is_first = exon_rank == 1,
    is_last  = exon_rank == max(exon_rank)
  ) %>%
  filter(is_first | is_last) %>%
  mutate(term_type = if_else(is_first, "five_prime", "three_prime")) %>%
  ungroup()

# Create the Snapping Map using the Mode
# We group by the 'Anchor' (the internal splice site) to find the consensus jitter site
snap_map <- terminal_exons %>%
  group_by(seqnames, strand, term_type, 
           anchor = if_else((strand == "+" & term_type == "five_prime") | 
                              (strand == "-" & term_type == "three_prime"), end, start)) %>%
  mutate(
    new_coord = {
      # Identify which side is jittery based on strand and type
      jitter_coords <- case_when(
        term_type == "five_prime" & strand == "+" ~ start,
        term_type == "five_prime" & strand == "-" ~ end,
        term_type == "three_prime" & strand == "+" ~ end,
        term_type == "three_prime" & strand == "-" ~ start
      )
      # Calculate Mode
      val_counts <- table(jitter_coords)
      modes <- as.numeric(names(val_counts[val_counts == max(val_counts)]))
      # Tie-breaker: pick the most upstream mode
      if(unique(strand) == "+") min(modes) else max(modes)
    }
  ) %>%
  ungroup() %>%
  dplyr::select(transcript_id, term_type, new_coord) %>%
  # Pivot so we have one row per transcript with both potential new coords
  pivot_wider(names_from = term_type, values_from = new_coord)

# Apply the Snapped Coordinates back to the GTF
gtf_step2 <- gtf_step2 %>%
  left_join(snap_map, by = "transcript_id") %>%
  group_by(transcript_id) %>%
  mutate(
    # Update Start
    start = case_when(
      type == "exon" & strand == "+" & start == min(start) & !is.na(five_prime) ~ five_prime,
      type == "exon" & strand == "-" & start == min(start) & !is.na(three_prime) ~ three_prime,
      TRUE ~ start
    ),
    # Update End
    end = case_when(
      type == "exon" & strand == "+" & end == max(end) & !is.na(three_prime) ~ three_prime,
      type == "exon" & strand == "-" & end == max(end) & !is.na(five_prime) ~ five_prime,
      TRUE ~ end
    )
  ) %>%
  # Sync the 'transcript' feature lines
  mutate(
    start = if_else(type == "transcript", min(start[type == "exon"]), start),
    end = if_else(type == "transcript", max(end[type == "exon"]), end)
  ) %>%
  ungroup() %>%
  dplyr::select(-five_prime, -three_prime)

# Clean up
rm(terminal_exons)

################################################################################
### (5) Correct 3' and 5' collapse + remove nested/truncated isoforms
################################################################################

# Collapse based on prefix/suffix junctions (5'/3' degradation)

# This version logs exact degradation type (5_prime_deg, 3_prime_deg, 5_3_prime_deg).
# Only compares isoforms with identical internal junctions.
# Detects 5′, 3′, or both terminal truncations.
# Keeps the longest isoform.
# Logs the removed isoforms for QC.

# Function to identify if one junction chain is a prefix or suffix of another
collapse_truncations <- function(gtf_df, log_removed_df = NULL) {
  
  # Prepare Junction Strings
  # We use the internal junctions to define the "backbone"
  tx_meta <- gtf_df %>%
    filter(type == "exon") %>%
    group_by(transcript_id, gene_id) %>%
    arrange(start) %>%
    summarise(
      # Junctions are the spaces between exons: end_i to start_i+1
      junc_string = if(n() > 1) paste0(end[-n()], "-", start[-1], collapse = ";") else "single-exon",
      n_exons = n(),
      .groups = "drop"
    )
  
  isoforms_to_remove <- data.frame(transcript_id = character(), reason = character())
  
  # Compare within each Gene
  for(g in unique(tx_meta$gene_id)) {
    gene_isoforms <- tx_meta %>% filter(gene_id == g) %>% arrange(desc(n_exons))
    
    if(nrow(gene_isoforms) < 2) next
    
    # We iterate through isoforms and see if they are "contained" in a longer one
    for(i in 1:nrow(gene_isoforms)) {
      query_id <- gene_isoforms$transcript_id[i]
      
      # If already marked for removal, skip
      if(query_id %in% isoforms_to_remove$transcript_id) next
      
      query_junc <- gene_isoforms$junc_string[i]
      
      # Compare against all other isoforms in the same gene
      for(j in 1:nrow(gene_isoforms)) {
        if(i == j) next # Don't compare to self
        
        target_id <- gene_isoforms$transcript_id[j]
        target_junc <- gene_isoforms$junc_string[j]
        
        # Logic: Is the query a suffix or prefix of the target?
        # A suffix starts with the same junction pattern but ends early (or vice versa)
        is_prefix <- str_starts(target_junc, fixed(query_junc))
        is_suffix <- str_ends(target_junc, fixed(query_junc))
        
        if((is_prefix | is_suffix) && nchar(query_junc) < nchar(target_junc)) {
          # The shorter one is the degradation product
          reason_str <- case_when(
            is_prefix & is_suffix ~ "internal_fragment",
            is_prefix ~ "3_prime_deg", # Target continues further 3', so query is 3' truncated
            is_suffix ~ "5_prime_deg"  # Target starts earlier 5', so query is 5' truncated
          )
          
          isoforms_to_remove <- rbind(isoforms_to_remove, 
                                      data.frame(transcript_id = query_id, reason = reason_str))
          break # Found a match, move to next query
        }
      }
    }
  }
  
  # Finalize
  gtf_filtered <- gtf_df %>% filter(!(transcript_id %in% isoforms_to_remove$transcript_id))
  
  if(!is.null(log_removed_df)) {
    log_removed_df <- rbind(log_removed_df, isoforms_to_remove)
  } else {
    log_removed_df <- isoforms_to_remove
  }
  
  return(list(gtf = gtf_filtered, log_removed = log_removed_df))
}

# Run the updated function
nested_results <- collapse_truncations(gtf_step2, log_removed)
gtf_step3 <- nested_results$gtf
log_removed <- nested_results$log_removed

# Filter ORF table to match surviving transcripts
orf_df_passed <- orf_df %>% 
  filter(transcript_id %in% gtf_step3$transcript_id)

################################################################################
### (x) Collapse ORF redundancy - DON'T DO
################################################################################

# # Join filtered classification data with expression data 
# orf_with_expr <- orf_df_passed %>%
#   full_join(fl_df, by = "transcript_id") %>%
#   mutate(fl_assoc = as.numeric(fl_assoc))  # ensure expression is numeric
# 
# # Break ties using transcript_id numeric suffix if needed
# orf_with_expr <- orf_with_expr %>%
#   mutate(
#     last_num = as.integer(str_extract(transcript_id, "(?<=\\.)\\d+$"))
#   )
# 
# # Separate isoforms into two groups:
# # 1) Those with ORF_seq non-NA (to collapse)
# to_collapse <- orf_with_expr %>%
#   filter(!is.na(orf_seq))
# # 2) Those with ORF_seq NA (keep as-is)
# to_keep_na <- orf_with_expr %>%
#   filter(is.na(orf_seq))
# # Add back gene_id
# to_keep_na$gene_id <- sub("^(([^.]+\\.[^.]+))\\..*$", "\\1", to_keep_na$transcript_id)
# 
# # For each gene + ORF_seq group, keep isoform with highest expression
# orf_collapse <- to_collapse %>%
#   group_by(gene_id, orf_seq) %>%
#   filter(fl_assoc == max(fl_assoc, na.rm = TRUE)) %>%
#   slice_min(order_by = last_num, with_ties = FALSE) %>%
#   ungroup()
# 
# # Combine collapsed non-NA isoforms with uncollapsed NA isoforms
# orf_collapse_na <- bind_rows(
#   orf_collapse,
#   to_keep_na
# )
# 
# # # Keep log of isoforms lost at this step
# # log_removed <- rbind(log_removed,
# #                      data.frame(
# #                        transcript_id = setdiff(orf_with_expr$transcript_id, orf_collapse_na$transcript_id),
# #                        reason = "ORF_redundancy",
# #                        stringsAsFactors = FALSE
# #                      ))
# # 
# # # Filter the GTF data frame to keep only rows matching those transcript IDs
# # gtf_step3 <- gtf_step2[gtf_step2$transcript_id %in% orf_collapse_na$transcript_id, ]
# 
# # Remove objects to clean R environment
# rm(to_collapse)
# rm(to_keep_na)
# rm(orf_with_expr)
# rm(orf_df)
# rm(orf_collapse)

################################################################################
### (x) Re-cluster within pb_ids to determine if there are multiple genes
################################################################################

# # CRITERIA FOR SEPARATING:
# # no shared intron junctions
# # less than 5% overlap
# 
# # Convert to data.table if not already
# if (!"data.table" %in% class(gtf_step2)) {
#   gtf_dt <- as.data.table(gtf_step2)
# } else {
#   gtf_dt <- copy(gtf_step2)
# }
# 
# all_pb_ids <- unique(gtf_dt$gene_id)
# names(all_pb_ids) <- all_pb_ids
# num_cores <- max(1, parallel::detectCores() - 1)
# cat("Using", num_cores, "cores for parallel processing\n")
# 
# # Define functions
# reciprocal_overlap <- function(start1, end1, start2, end2) {
#   overlap <- max(0, min(end1, end2) - max(start1, start2) + 1)
#   len1 <- end1 - start1 + 1
#   len2 <- end2 - start2 + 1
#   return(overlap / min(len1, len2))
# }
# 
# get_junctions <- function(exons) {
#   # Ensure exons are sorted by start
#   exons <- exons[order(exons$start), ]
#   # If fewer than 2 exons, no junctions
#   if (nrow(exons) < 2) return(character(0))
#   # Junctions: end of exon i and start of exon i+1
#   junctions <- paste0(exons$end[-nrow(exons)], "_", exons$start[-1])
#   return(junctions)
# }
# 
# cluster_one_gene <- function(pb_id) {
#   df_group <- gtf_dt[gene_id == pb_id]
#   transcripts <- unique(df_group$transcript_id)
#   
#   if (length(transcripts) == 1) {
#     df_group[, new_gene_id := paste0(pb_id, "_1")]
#     return(list(data = df_group, count = 1))
#   }
#   
#   # Precompute bounds and junctions per transcript
#   transcript_info <- lapply(transcripts, function(tx) {
#     exons <- df_group[transcript_id == tx & type == "exon"]
#     list(
#       tx = tx,
#       start = min(exons$start),
#       end = max(exons$end),
#       junctions = get_junctions(exons)
#     )
#   })
#   
#   edges_list <- list()
#   n <- length(transcripts)
#   
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       t1 <- transcript_info[[i]]
#       t2 <- transcript_info[[j]]
#       
#       shared <- length(intersect(t1$junctions, t2$junctions)) > 0
#       r_overlap <- reciprocal_overlap(t1$start, t1$end, t2$start, t2$end)
#       
#       if (shared || r_overlap > 0.05) {
#         edges_list <- append(edges_list, list(c(t1$tx, t2$tx)))
#       }
#     }
#   }
#   
#   # Always create graph with all transcripts as vertices
#   g <- igraph::make_empty_graph(n = length(transcripts), directed = FALSE)
#   g <- igraph::set_vertex_attr(g, "name", value = transcripts)
#   
#   if (length(edges_list) > 0) {
#     edges_mat <- do.call(rbind, edges_list)
#     edges_vec <- as.vector(t(edges_mat))
#     g <- igraph::add_edges(g, edges_vec)
#   }
#   
#   comps <- igraph::components(g)
#   
#   new_ids <- paste0(pb_id, "_", comps$membership[igraph::V(g)$name])
#   tx_to_new_gene <- setNames(new_ids, igraph::V(g)$name)
#   
#   df_group[, new_gene_id := tx_to_new_gene[transcript_id]]
#   
#   return(list(data = df_group, count = length(unique(new_ids))))
# }
# 
# # Run in parallel
# results <- mclapply(all_pb_ids, 
#                     cluster_one_gene, 
#                     mc.cores = num_cores)
# 
# # Combine results
# reclustered_rows <- rbindlist(lapply(results, `[[`, "data"))
# cluster_counts <- sapply(results, `[[`, "count")
# 
# cat("Done reclustering", length(all_pb_ids), "genes\n")
# cat("Mean clusters per gene:", mean(cluster_counts), "\n")
# cat("Median clusters per gene:", median(cluster_counts), "\n")
# 
# # Combine into final data.frame
# gtf_reclustered <- bind_rows(reclustered_rows)
# 
# # Report summary
# summary_df <- tibble(
#   gene_id = names(cluster_counts),
#   new_gene_count = unlist(cluster_counts)
# )
# 
# #### Sort gtf
# gtf_sorted <- gtf_reclustered %>%
#   arrange(seqnames, start, end, strand)
# 
# unique_genes <- unique(gtf_sorted$new_gene_id)
# new_gene_ids <- paste0("G.", seq_along(unique_genes))
# 
# gene_id_map <- setNames(new_gene_ids, unique_genes)
# 
# # Replace new_gene_id with new G.# IDs
# gtf_sorted <- gtf_sorted %>%
#   mutate(new_new_gene_id = gene_id_map[new_gene_id])
# 
# # Map old transcript_ids to new transcript numbers within each new_new_gene_ids
# transcript_map <- gtf_sorted %>%
#   dplyr::select(new_new_gene_id, transcript_id) %>%
#   distinct() %>%
#   group_by(new_new_gene_id) %>%
#   arrange(transcript_id) %>%  # or arrange by something else if needed
#   mutate(transcript_num = row_number()) %>%
#   ungroup() %>%
#   mutate(new_transcript_id = paste0(new_new_gene_id, ".", transcript_num))
# 
# # Join new transcript IDs back to the main dataframe
# gtf_sorted <- gtf_sorted %>%
#   left_join(transcript_map %>% 
#               dplyr::select(transcript_id, new_transcript_id), by = "transcript_id") 
# 
# # Create summary dataframe
# gene_reassignment <- gtf_sorted %>%
#   dplyr::select(pb_id = gene_id,
#          transcript_id,
#          new_gene_id = new_new_gene_id,
#          new_transcript_id) %>%
#   distinct()
# 
# # Save to file (e.g., CSV)
# write.csv(gene_reassignment,
#           here::here("data",
#                      "isoseq",
#                      "gene_reassignment_map.csv"),
#           row.names = FALSE)
# 
# # Clean up gtf
# gtf_final <- gtf_sorted %>%
#   mutate(gene_id = new_new_gene_id,
#          transcript_id = new_transcript_id) %>%
#   dplyr::select(-new_new_gene_id, -new_transcript_id, -new_gene_id)
# 
# # Remove objects to tidy R environment
# rm(all_pb_ids)
# rm(num_cores)
# rm(cluster_one_gene)
# rm(get_junctions)
# rm(gtf_sorted)

################################################################################
### (6) Switch gene names to reference names using SQANTI output
################################################################################

# Load sqanti results with reference name
asso_gene <- read.delim2(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "12_sqanti_filter",
                                 "SQANTI_filtered_RulesFilter_result_classification.txt")) %>%
  filter(filter_result == "Isoform") %>%
  dplyr::select(associated_gene, isoform)

glimpse(gtf_step3)
glimpse(asso_gene)

# Replace gtf_step3$gene_id with asso_gene$associated_gene using transcript_id == isoform
gtf_step3 <- gtf_step3 %>%
  # Bring in the associated_gene from SQANTI
  left_join(asso_gene, by = c("transcript_id" = "isoform")) %>%
  # Create new_gene_id (fallback to old gene_id if no SQANTI match)
  mutate(new_gene_id = coalesce(associated_gene, gene_id)) %>%
  # Create new_transcript_id with version suffix
  # Grouping ensures the numbering starts over for each gene
  group_by(new_gene_id) %>%
  mutate(
    # match() + unique() creates a sequence (1, 2, 3...) for each unique transcript within the gene
    version_num = match(transcript_id, unique(transcript_id)),
    new_transcript_id = paste0(new_gene_id, ".", str_pad(version_num, 1))
  ) %>%
  ungroup() %>%
  # Clean up intermediate columns
  dplyr::select(-associated_gene, -version_num)

# Find all genes/transcripts with '_' in them, representing fusion genes
# When I look at these in IGV they don't look like fusions, but there are enough of them
# n=675 transcript that I don't want to manually check them, so I'm just going to remove them
underscore_ids <- gtf_step3 %>%
  filter(grepl("_", new_gene_id)) %>%
  distinct(new_gene_id)

transcripts_to_remove <- gtf_step3 %>%
  filter(grepl("_", new_gene_id)) %>%
  pull(new_transcript_id) %>%
  unique()

cat("Found", nrow(underscore_ids), "gene IDs with cryptic fusions")
cat("Found", length(transcripts_to_remove), "transcript IDs with cryptic fusions")

# Filter the original dataframe to REMOVE those rows
gtf_step3 <- gtf_step3 %>%
  filter(!grepl("_", new_gene_id))

# Log them safely
if (length(transcripts_to_remove) > 0) {
  log_removed <- rbind(
    log_removed,
    data.frame(
      transcript_id = transcripts_to_remove,
      reason = rep("non_coding_gene", length(transcripts_to_remove)),
      stringsAsFactors = FALSE
    )
  )
}

################################################################################
### (7) Remove non-coding genes
################################################################################

# Update your ORF data to use the new naming convention
orf_joined <- orf_df_passed %>%
  # We select only the 'new' ID columns from the GTF to bring into the ORF data
  right_join(
    gtf_step3 %>% 
      dplyr::select(gene_id, transcript_id, new_gene_id, new_transcript_id) %>% 
      distinct(), 
    by = c("gene_id", "transcript_id") # Join on the original PacBio IDs
  ) %>%
  # Now overwrite the standard columns with your 'new' names for the audit
  mutate(
    old_transcript_id = transcript_id,
    old_gene_id = gene_id,
    transcript_id = new_transcript_id,
    gene_id = new_gene_id
  )

glimpse(orf_joined)

# Summarize the 'Best' ORF status for every gene
gene_audit <- orf_joined %>%
  group_by(new_gene_id) %>%
  summarise(
    # This will now be FALSE for your 1,296 missing genes
    has_any_orf = any(!is.na(orf_seq)), 
    has_complete = any(orf_type == "complete", na.rm = TRUE),
    only_incomplete = has_any_orf & !has_complete,
    .groups = "drop"
  )

glimpse(gene_audit)

# Get the count of genes that fall into the 'Incomplete Only' trap
incomplete_only_count <- sum(gene_audit$only_incomplete)
print(paste("Genes with ONLY incomplete ORFs:", incomplete_only_count))

# Identify the gene_ids to remove
genes_to_remove <- gene_audit %>%
  filter(has_any_orf == FALSE | only_incomplete == TRUE) %>%
  pull(new_gene_id)

# Get transcript IDs belonging to genes marked for removal
transcripts_to_remove <- gtf_step3 %>% 
  filter(new_gene_id %in% genes_to_remove) %>%
  pull(transcript_id) %>% 
  unique()

# Final Filter of the GTF
gtf_final <- gtf_step3 %>%
  filter(!new_gene_id %in% genes_to_remove)

# Log them safely
if (length(transcripts_to_remove) > 0) {
  log_removed <- rbind(
    log_removed,
    data.frame(
      transcript_id = transcripts_to_remove,
      reason = rep("non_coding_gene", length(transcripts_to_remove)),
      stringsAsFactors = FALSE
    )
  )
}

cat("Removed", length(genes_to_remove), "genes for being entirely non-coding.\n")

# Remove objects to tidy R environment
rm(genes_to_remove)

################################################################################
### (8) Remove non-autosomal genes
################################################################################

# Identify non-autosomal genes
non_autosomal_genes <- gtf_final %>%
  filter(
    grepl("^NW", as.character(seqnames)) | # unplaced
      seqnames == "NC_024342.1" | # sex chromosome
      seqnames == "NC_024238.1" # remove mitochondria# Specific sex chromosome
  ) %>%
  pull(gene_id) %>%
  unique()

# Filter to keep only autosomal genes
gtf_final2 <- gtf_final %>%
  filter(!(gene_id %in% non_autosomal_genes))

# Log them safely
if (length(transcripts_to_remove) > 0) {
  log_removed <- rbind(
    log_removed,
    data.frame(
      transcript_id = transcripts_to_remove,
      reason = rep("non_coding_gene", length(transcripts_to_remove)),
      stringsAsFactors = FALSE
    )
  )
}

cat("Removed", length(non_autosomal_genes), "genes for being non-autosomal.\n")

# Remove objects to tidy R environment
rm(non_autosomal_genes)

################################################################################
### (9) Save results
################################################################################

# Create valid GTF format from your df
gtf_formatted <- data.frame(
  seqname   = as.character(gtf_final2$seqnames),  # convert factor to character
  source    = as.character(gtf_final2$source),    # convert factor to character
  feature   = as.character(gtf_final2$type),      # convert factor to character
  start     = gtf_final2$start,
  end       = gtf_final2$end,
  score     = ifelse(is.na(gtf_final2$score), ".", gtf_final2$score),
  strand    = as.character(gtf_final2$strand),    # convert factor to character
  frame     = ifelse(is.na(gtf_final2$phase), ".", gtf_final2$phase),
  attribute = paste0(
    'gene_id "', gtf_final2$new_gene_id, '"; ',
    'transcript_id "', gtf_final2$new_transcript_id, '";'
  ),
  stringsAsFactors = FALSE
)

# Save new gtf!
write.table(
  gtf_formatted,
  file = here::here("data",
                    "GTF",
                    "03_nonredundant_transcriptome.gtf"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Keep a copy in this step's folder too
file.copy(here::here("data", "GTF", "03_nonredundant_transcriptome.gtf"),
          here::here("data", "01_build_longread_transcriptome", "14_remove_redundancy", "03_nonredundant_transcriptome.gtf"),
          overwrite = TRUE)

# Print a summary
cat("Non-redundant gtf save with", length(unique(gtf_final2$new_gene_id)), "unique genes.\n")
cat("Non-redundant gtf save with", length(unique(gtf_final2$new_transcript_id)), "unique transcripts.\n")

################################################################################
### (10) Filter fasta to match
################################################################################

fasta <- readDNAStringSet(here::here("data",
                                     "01_build_longread_transcriptome",
                                     "12_sqanti_filter",
                                     "SQANTI_filtered.filtered.fasta"))

# Remove sequences that are in log_removed
fasta_filtered <- fasta[!(names(fasta) %in% log_removed$transcript_id)]

# Create a named vector to map old names to new ones
name_map_df <- gtf_final2 %>%
  dplyr::select(transcript_id, new_transcript_id) %>%
  distinct()

name_map <- setNames(name_map_df$new_transcript_id, name_map_df$transcript_id)

# Only keep sequences that exist in our final filtered GTF 
# (This catches anything dropped in the ORF filter)
fasta_final <- fasta_filtered[names(fasta_filtered) %in% names(name_map)]

# Apply the renaming
names(fasta_final) <- name_map[names(fasta_final)]

# Save
writeXStringSet(fasta_final, filepath = here::here("data",
                                                      "01_build_longread_transcriptome",
                                                      "14_remove_redundancy",
                                                      "nonredundant_transcripts.fasta"))

# Save map 'rosetta stone'
# Create a clean data frame mapping old IDs to new IDs
id_mapping_table <- gtf_final2 %>%
  dplyr::select(gene_id, transcript_id, new_gene_id, new_transcript_id) %>%
  distinct()

# Save as a CSV for easy opening in Excel or R later
write.csv(
  id_mapping_table, 
  file = here::here("data",
                    "01_build_longread_transcriptome",
                    "14_remove_redundancy",
                    "rosetta_stone.csv"),
  row.names = FALSE
)

# Print a summary
cat("Mapping file saved with", length(unique(id_mapping_table$new_gene_id)), "unique genes.\n")
cat("Mapping file saved with", nrow(id_mapping_table), "unique transcripts.\n")

