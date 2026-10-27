### Filter stringtie merged.gtf before kallisto

### FILTER NOTICE: remove any transcript with unsupported STAR junctions
### FILTER NOTICE: remove non-autosomal
### FILTER NOTICE: remove novel genes
### FILTER NOTICE: remove all mono-exonic isoforms
### FILTER NOTICE: remove non-coding genes
### FILTER NOTICE: collapse isoforms the share all internal junctions

library(rtracklayer)
library(dplyr)
library(stringr)

################################################################################
### (1) Load
################################################################################

# Load the GTF file
stringtie_gtf <- rtracklayer::import(here::here("data",
                                                "02_build_shortread_transcriptome",
                                                "01_stringtie",
                                                "merge",
                                                "merged.gtf"))

stringtie_df <- data.frame(stringtie_gtf)

# Print a summary
cat("Stringtie GTF contains", length(unique(stringtie_df$ref_gene_id)), "unique ref genes.\n")
cat("Stringtie GTF contains", length(unique(stringtie_df$gene_id)), "unique genes.\n")
cat("Stringtie GTF contains", length(unique(stringtie_df$transcript_id)), "unique transcripts.\n")

# Remove unstranded rows from the dataframe
# Identify IDs of transcripts that have at least one "*" row
unstranded_ids <- stringtie_df %>%
  filter(strand == "*") %>%
  pull(transcript_id) %>%
  unique()

# Remove those entire transcripts from the dataframe
stringtie_df_filtered <- stringtie_df %>%
  filter(!(transcript_id %in% unstranded_ids))

# Print a clean summary
cat("Identified", length(unstranded_ids), "transcripts with unstranded features.\n")
cat("Removed a total of", nrow(stringtie_df) - nrow(stringtie_df_filtered), "rows (transcripts + exons).\n")

cat("Stringtie GTF contains", length(unique(stringtie_df_filtered$ref_gene_id)), "unique ref genes.\n")
cat("Stringtie GTF contains", length(unique(stringtie_df_filtered$gene_id)), "unique genes.\n")
cat("Stringtie GTF contains", length(unique(stringtie_df_filtered$transcript_id)), "unique transcripts.\n")

sj_supported <- read.csv(here::here("data",
                                    "01_build_longread_transcriptome",
                                    "14_remove_redundancy",
                                    "junction_whitelist.csv"))

################################################################################
### (2) Filter by supported junctions
################################################################################

# Identify all multi-exon transcripts and their junctions
stringtie_junctions <- stringtie_df_filtered %>%
  filter(type == "exon") %>%
  group_by(transcript_id) %>%
  filter(n() > 1) %>% # Only transcripts with >1 exon have junctions
  arrange(start) %>%
  mutate(
    intron_start = end + 1,
    intron_end = lead(start) - 1
  ) %>%
  filter(!is.na(intron_end)) %>% # Remove the last row of each group
  dplyr::select(transcript_id, seqnames, intron_start, intron_end, strand) %>%
  ungroup()

# Create the whitelist match string (Chr-Start-End-Strand)
whitelist_keys <- sj_supported %>%
  mutate(match_key = paste(chr, intron_start, intron_end, strand_label, sep = "_")) %>%
  pull(match_key)

# Check every junction against the whitelist
stringtie_junctions <- stringtie_junctions %>%
  mutate(
    match_key = paste(seqnames, intron_start, intron_end, strand, sep = "_"),
    is_supported = match_key %in% whitelist_keys
  )

# Identify Transcripts to Remove
# A transcript is "bad" if ANY of its junctions are NOT supported
failed_tx_ids <- stringtie_junctions %>%
  group_by(transcript_id) %>%
  summarise(all_passed = all(is_supported)) %>%
  filter(!all_passed) %>%
  pull(transcript_id)

# Final Filter
final_stringtie_df <- stringtie_df_filtered %>%
  filter(!(transcript_id %in% failed_tx_ids))

# Summary Stats
cat("Transcripts removed (unstranded):", length(unstranded_ids), "\n")
cat("Transcripts removed (unsupported junctions):", length(failed_tx_ids), "\n")
cat("Transcripts remaining:", length(unique(final_stringtie_df$transcript_id)), "\n")
cat("Genes remaining:", length(unique(final_stringtie_df$gene_id)), "\n")
cat("Ref genes remaining:", length(unique(final_stringtie_df$ref_gene_id)), "\n")

################################################################################
### (3) Remove non-autosomal genes
################################################################################

# Define Autosomal Chromosomes (based on your earlier grep filter)
# This removes NC_024342.1 and any unplaced scaffolds starting with "NW"
autosomal_df <- final_stringtie_df %>%
  filter(seqnames != "NC_024342.1") %>% # remove sex chromosome
  filter(!str_detect(seqnames, "^NW")) %>% # remove unplaced scaffolds
  filter(seqnames != "NC_024238.1") %>% # remove mitochondria
  filter(!str_detect(ref_gene_id, "^unassigned")) # remove novel genes

cat("Transcripts remaining:", length(unique(autosomal_df$transcript_id)), "\n")
cat("Genes remaining:", length(unique(autosomal_df$gene_id)), "\n")
cat("Ref genes remaining:", length(unique(autosomal_df$ref_gene_id)), "\n")

################################################################################
### (4) Remove putative fusions/ambiguous mapping
################################################################################

# Identify fusion genes (gene_id mapping to multiple ref_gene_id)
fusion_gene_ids <- autosomal_df %>%
  group_by(gene_id) %>%
  summarise(n_ref_gene = n_distinct(ref_gene_id), .groups = "drop") %>%
  filter(n_ref_gene > 2) %>% # 2 so it can include NAs
  pull(gene_id)

cat("Number of putative fusion genes to remove:", length(fusion_gene_ids), "\n")

# Remove these from your dataframe
no_fusions <- autosomal_df %>%
  filter(!(gene_id %in% fusion_gene_ids))

# Check counts
cat("Transcripts remaining:", length(unique(no_fusions$transcript_id)), "\n")
cat("Genes remaining:", length(unique(no_fusions$gene_id)), "\n")
cat("Ref genes remaining:", length(unique(no_fusions$ref_gene_id)), "\n")

# Assign ref_gene_id to gene_id
no_fusions$gene_id <- no_fusions$ref_gene_id

# Check counts
cat("Transcripts remaining:", length(unique(no_fusions$transcript_id)), "\n")
cat("Genes remaining:", length(unique(no_fusions$gene_id)), "\n")
cat("Ref genes remaining:", length(unique(no_fusions$ref_gene_id)), "\n")

################################################################################
### (5) Assign ref_gene_id and remove novel genes - no longer needed
################################################################################
# 
# ### Propogate ref_gene_id to transcripts that don't have it
# # but match gene_id and overlap
# # We group by the StringTie gene_id (MSTRG.#) to find "siblings"
# no_unassigned <- no_fusions %>%
#   group_by(gene_id) %>%
#   mutate(
#     # Find the first 'real' ID in the group that doesn't start with "unassigned"
#     true_ref_id = ref_gene_id[!str_detect(ref_gene_id, "^unassigned")][1],
#     # If a true ID exists for this MSTRG cluster, assign it to everyone
#     ref_gene_id = if_else(!is.na(true_ref_id), true_ref_id, ref_gene_id),
#   ) %>%
#   select(-true_ref_id) %>%
#   ungroup()
# 
# # Filter out Novel/Unassigned Genes
# # We keep only records with a valid ref_gene_id and 
# # exclude StringTie's "unassigned_gene" placeholders.
# known_genes_df <- no_unassigned %>%
#   filter(!is.na(ref_gene_id),
#          !str_detect(ref_gene_id, "^unassigned"))
# 
# cat("Transcripts remaining:", length(unique(known_genes_df$transcript_id)), "\n")
# cat("Genes remaining:", length(unique(known_genes_df$gene_id)), "\n")
# cat("Ref genes remaining:", length(unique(known_genes_df$ref_gene_id)), "\n")

################################################################################
### (6) Final mono-exonic genes
################################################################################

# Final Single-Exon Filter (Kallisto/TransDecoder safety)
# A transcript needs a 'transcript' row PLUS at least one 'exon' row.
final_filtered_df <- no_fusions %>%
  group_by(transcript_id) %>%
  filter(any(type == "exon")) %>% 
  filter(n() > 2) %>%             
  ungroup()

cat("Transcripts remaining:", length(unique(final_filtered_df$transcript_id)), "\n")
cat("Genes remaining:", length(unique(final_filtered_df$gene_id)), "\n")
cat("Ref genes remaining:", length(unique(final_filtered_df$ref_gene_id)), "\n")

# Fix gene_id for ALL rows (Exons + Transcripts)
# We map the transcript_id back to our dataframe to get the correct ref_gene_id
# This ensures exons inherit the gene name, preventing TransDecoder crashes.
id_map <- final_filtered_df %>%
  dplyr::select(transcript_id, ref_gene_id) %>%
  distinct()

################################################################################
### (7) Remove non-coding genes
################################################################################

# Load td2 results
td2 <- read.csv(here::here("data",
                           "02_build_shortread_transcriptome",
                           "02_transdecoder",
                           "parsed_results.csv"))

# Remove all genes with no coding transcripts
# Summarize the 'Best' ORF status for every gene
gene_audit <- final_filtered_df %>%
  left_join(td2, by = join_by(transcript_id == transcript, ref_gene_id == gene)) %>%
  group_by(ref_gene_id) %>%
  summarise(
    has_any_orf = any(!is.na(orf_seq)),
    # Check if the gene has at least one complete ORF
    has_complete = any(orf_type == "complete", na.rm = TRUE),
    # Identify genes that only have incomplete ORFs
    only_incomplete = has_any_orf & !has_complete,
    .groups = "drop"
  )

# Identify genes that have at least one ORF (coding) and passed the audit
valid_genes <- gene_audit %>%
  filter(has_any_orf) %>%          # keep genes with any ORF
  pull(ref_gene_id)                # get vector of gene IDs

# Filter df
final_filtered_df2 <- final_filtered_df %>%
  filter(ref_gene_id %in% valid_genes)

# Check
cat("Filtered df contains", length(unique(final_filtered_df2$gene_id)), "unique genes.\n")
cat("Filtered df contains", length(unique(final_filtered_df2$ref_gene_id)), "unique ref genes.\n")
cat("Filtered df contains", length(unique(final_filtered_df2$transcript_id)), "unique transcripts.\n")

################################################################################
### (8) Collapse redundancy due to potential 5' 3' deg
################################################################################

# Collapse all genes that share internal junctions but are just degrades on one side

# Function to identify if one junction chain is a prefix or suffix of another
# from long-reads
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

# Run the collapse function on your short-read dataframe
short_read_results <- collapse_truncations(final_filtered_df2)

# Extract the clean GTF and the log
final_short_gtf_collapsed <- short_read_results$gtf
short_log_removed <- short_read_results$log_removed

# Summary Stats for the Short-Read Cleanup
cat("Genes after collapse:", length(unique(final_short_gtf_collapsed$gene_id)), "\n")
cat("Transcripts after collapse:", length(unique(final_short_gtf_collapsed$transcript_id)), "\n")
cat("Redundant isoforms removed:", nrow(short_log_removed), "\n")

################################################################################
### (9) Save
################################################################################

glimpse(final_short_gtf_collapsed)

# Change gtf so that gene_id is ref_gene_id
mcols(stringtie_gtf)$gene_id <- id_map$ref_gene_id[match(stringtie_gtf$transcript_id, id_map$transcript_id)]

# Vector of valid genes and transcripts
final_valid_transcripts <- unique(final_short_gtf_collapsed$transcript_id)
valid_genes <- unique(final_short_gtf_collapsed$ref_gene_id)

# Filter the GTF
filtered_gtf_object <- stringtie_gtf[mcols(stringtie_gtf)$transcript_id %in% final_valid_transcripts]

# Check counts
cat("Transcripts remaining:", length(unique(mcols(filtered_gtf_object)$transcript_id)), "\n")
cat("Genes remaining:", length(unique(mcols(filtered_gtf_object)$gene_id)), "\n")

# Export
rtracklayer::export(filtered_gtf_object, 
                    here::here("data", 
                               "02_build_shortread_transcriptome", 
                               "03_filter_stringtie",
                               "merged_filtered.gtf"), 
                    format = "gtf")


# FINAL FILTERING RESULTS:
# step	gene	ref_gene	transcript
# stringtie	29399	26259	82395
# remove unstranded	24466	26259	77462
# remove unsupported	21587	21394	55374
# remove non-autosomal	18036	19330	32924
# remove fusions	17887	18832	32098
# assign ref_ids	18832	18832	32098
# remove mono-exon	16960	16960	29901
# remove non-coding	14269	14269	26358
# remove degredation	14269	14269	25284