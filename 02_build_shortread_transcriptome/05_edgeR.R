### Use edgeR to get scaled and filtered count matrices from SRP

# Load data
# Use edgeR to load counts, scale and filter

library(tidyr)
library(ggplot2)
library(dplyr)
library(edgeR)
library(rtracklayer)

################################################################################
################################################################################
### SHORT READS
################################################################################
################################################################################

################################################################################
### (1) Count matrices
################################################################################

# Get kallisto output paths

# edgeR wants the directory where the abundance.h5 files are kept as an input

# # Load gtf
gtf <- import(here::here("data",
                          "GTF",
                          "05_shortreadpipe.gtf"))
gtf_df <- as.data.frame(gtf)

# Check counts
cat("Genes:", length(unique(gtf_df$gene_id)), "\n")
cat("Transcripts:", length(unique(gtf_df$transcript_id)), "\n")

# List parent directory
base_dir <- here::here("data",
                       "02_build_shortread_transcriptome",
                       "04_kallisto",
                       "kallisto_short")

# Get all sample directory names
sample_dirs <- list.dirs(base_dir, 
                         recursive = FALSE)

# Build the full paths to abundance_1 folders
abundance_dirs <- file.path(sample_dirs, 
                            "quant_unfiltered", 
                            "abundance_1")

# Load transcript-level reads with catchKallisto

# Read in the data using catchKallisto
y <- catchKallisto(abundance_dirs)

 ### FILTER NOTICE: Remove all non-autosomal chromosomes (this is just a double check)
 # Remove non-autosomal genes (Chr 21 and unplaced scaffolds)
 gtf_df_filtered <- gtf_df[!grepl("^NW", as.character(gtf_df$seqnames)), ] %>%
   filter(seqnames != "NC_024342.1") %>% # remove transcripts on 12
   filter(type == "transcript") # only keep transcript rows
 
 # Also load metadata for batch effects
 # Contains: cell, sample_id, sex, tissue, selection, read depth (from flnc.bam)
 metadata <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "metadata.csv"))
 
# # Filter GTF for transcripts without gene_id
# gtf_df_filtered <- gtf_df_filtered %>%
#   filter(type == "transcript" & !is.na(gene_id))
# # 
# # Load td2 results
# td2 <- read.csv(here::here("data",
#                            "isoseq",
#                            "14_shortreadpipe",
#                            "04_TD2",
#                            "parsed_results.csv"))
# 
# # Remove all genes with no coding transcripts
# # Summarize the 'Best' ORF status for every gene
# gene_audit <- gtf_df_filtered %>%
#   left_join(td2, by = join_by(transcript_id == transcript, gene_id == gene)) %>%
#   group_by(gene_id) %>%
#   summarise(
#     has_any_orf = any(!is.na(orf_seq)),
#     # Check if the gene has at least one complete ORF
#     has_complete = any(orf_type == "complete", na.rm = TRUE),
#     # Identify genes that only have incomplete ORFs
#     only_incomplete = has_any_orf & !has_complete,
#     .groups = "drop"
#   )
# 
# # Get the count of genes that fall into the 'Incomplete Only' trap
# incomplete_only_count <- sum(gene_audit$only_incomplete)
# print(paste("Genes with ONLY incomplete ORFs:", incomplete_only_count))
# 
# # Identify the gene_ids to remove
# genes_to_remove <- gene_audit %>%
#   filter(has_any_orf == FALSE | only_incomplete == TRUE) %>%
#   pull(gene_id)
# 
# ### FILTER NOTICE: Remove non-coding genes
# # Final Filter of the GTF
# gtf_df_filtered <- gtf_df_filtered %>%
#   filter(!gene_id %in% genes_to_remove)# Final Filter of the GTF

# Filter
y[c("counts", "annotation")] <- lapply(
  y[c("counts", "annotation")],
  function(x) x[rownames(x) %in% unique(gtf_df_filtered$transcript_id), , drop = FALSE]
)

# Create metadata
# Your vector of all sample names
all_samples <- c("RFG1", "RFG2", "RFG3", "RFG4", "RFG5", "RFG6",
                 "RFT1", "RFT2", "RFT3", "RFT4", "RFT5", "RFT6",
                 "RMG1", "RMG2", "RMG3", "RMG4", "RMG5", "RMG6",
                 "RMT1", "RMT2", "RMT3", "RMT4", "RMT5", "RMT6")

# Combine into a metadata data frame
sample_info <- data.frame(
  sample = all_samples,
  sex = sex <- factor(ifelse(grepl("F", all_samples), "Female", "Male")),
  tissue = tissue <- factor(ifelse(grepl("T", all_samples), "Tail", "Gonad")),   # <-- FIXED here
  treatment = factor(paste0(sex, "_", tissue))
)

# Add batch effects from metadata (cell) and selection column
sample_info <- sample_info %>%
  left_join(metadata %>% 
              dplyr::select(sample_id, selection, cell),
            by = c("sample" = "sample_id"))

# Create transcript-level DGE objects

# Use divided count strategy to control for Overdispersion at transcript level
dge <- DGEList(counts = y$counts / y$annotation$Overdispersion,
               genes = y$annotation,
               samples = sample_info,
               group = sample_info$treatment)

### Add sample info to dge object (fixing directory names)
rownames(sample_info) <- sample_info$sample
# Replace column names with the part after flnc_ and before the next "\\"
colnames(dge) <- sub(".*/(R[FM][GT][0-9])/.*", "\\1", colnames(dge))
# Check the new column names
print(colnames(dge))
all(rownames(sample_info) == colnames(dge))  # should return TRUE
dge$samples <- cbind(sample_info[colnames(dge), ], 
                     dge$samples[, c("lib.size", "norm.factors")])

# Collapse to gene-level

# Create transcript-to-gene mapping from gtf_df
t2g <- unique(gtf_df[, c("transcript_id", "gene_id")])

# Add gene id and transcript id to dge
dge$genes$TranscriptID <- rownames(dge)
dge$genes$GeneID <- t2g$gene_id[match(rownames(dge), t2g$transcript_id)]

# Aggregate counts to gene-level
gene_counts <- rowsum(dge$counts, group = dge$genes$GeneID)

# Create gene-level annotation by taking the first transcript info per gene
gene_annotation <- dge$genes[!duplicated(dge$genes$GeneID), ]
rownames(gene_annotation) <- gene_annotation$GeneID

# Create dge for gene level
dge_gene <- DGEList(counts = gene_counts,
                    genes = gene_annotation[rownames(gene_counts), ],
                    samples = dge$samples,
                    group = dge$samples$treatment)

# Remove columns we no longer need (they are)
dge_gene$genes$TranscriptID <- NULL
dge_gene$samples$norm.factors.1 <- NULL
dge_gene$samples$lib.size.1 <- NULL

# Look at gene-level cpm distribution

# Calculate log2 CPM with prior.count = 1 (do this before normalizing)
logCPM <- cpm(dge_gene, log=TRUE, prior.count=0.1)

# Convert to data frame and add gene IDs
df <- as.data.frame(logCPM)
df$GeneID <- rownames(logCPM)

# Gather into long format: columns = GeneID, Sample, logCPM
df_long <- df %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "logCPM")

# Plot
ggplot(df_long, aes(x = logCPM, color = Sample)) +
  geom_density() +
  labs(title = "Distribution of log2 CPM values per sample",
       x = "Log2 CPM",
       y = "Density") +
  theme_minimal()

# Save 
ggsave(here::here("figures", 
                  "05_SRP_histogram_genecpm_short.png"), 
       width = 12, 
       height = 10, 
       dpi = 300)

### USE 1 AS A CPM CUT-OFF

# Filter dge_gene_gonad across 4 treatments and save

# Find genes to keep across all 4 treatments (used for other analyses)
# Criteria: CPM > 1 in at least 3/6 reps in at least 1 treatment
# DO THIS BEFORE NORMALIZING

# Logical matrix to store whether each gene passes in any treatment group
keep <- rep(FALSE, nrow(cpm(dge_gene)))  # start with all FALSE

# Calculate Total CPM once (average across all samples)
# This ensures the gene has a baseline presence across the whole experiment
total_cpm <- rowSums(cpm(dge_gene))

for (grp in unique(dge_gene$samples$treatment)) {
  idx <- which(dge_gene$samples$treatment == grp)
  # Check for CPM > 1 in at least 3 samples within the group
  # (Note: You mentioned TPM > 5 in comments but used > 1 in code; I kept > 1)
  keep_grp <- rowSums(cpm(dge_gene)[, idx] > 1) >= 3
  # Update global keep: 
  # Must pass the group rule AND the total CPM > 1 rule
  keep <- keep | (keep_grp & total_cpm > 0)
}

table(keep)
#FALSE  TRUE 
#324 13945 

# Subset and recalculate library sizes
dge_gene_4 <- dge_gene[keep, , keep.lib.sizes = FALSE]

# Remove transcripts with zero counts across all samples
nonzero_rows <- rowSums(dge_gene_4$counts) > 0
dge_gene_4 <- dge_gene_4[nonzero_rows, , keep.lib.sizes = FALSE]

# Normalize
dge_gene_4_norm <- calcNormFactors(dge_gene_4, method = "TMM")

# Export count matrix and TMM-CPM matrix
gene_matrix_4 <- as.data.frame(dge_gene_4_norm$counts)
gene_matrix_TPM4 <- as.data.frame(cpm(dge_gene_4_norm, normalized.lib.sizes = TRUE), log = F)

# Add gene column
gene_matrix_4$gene <- rownames(gene_matrix_4)
gene_matrix_CPM4$gene <- rownames(gene_matrix_CPM4)

# # Save
# write.csv(gene_matrix_4, 
#           here::here("data",
#                      "isoseq",
#                      "10_edgeR_filterexp",
#                      "01_all_gene_countmatrix_short.csv"))
# write.csv(gene_matrix_CPM4,
#           here::here("data",
#                      "isoseq",
#                      "10_edgeR_filterexp",
#                      "01_all_gene_TMMTPMmatrix_short.csv"))

# Subset DGE into each pairwise comparison in new objects and filter

# Helper function: Filter by treatment groups within a DGE subset
filter_by_treatment <- function(dge_obj, TPM_cutoff = 1, min_samples = 3, total_cutoff = 1) {
  TPM_matrix <- cpm(dge_obj)
  treatment_groups <- unique(dge_obj$samples$treatment)
  # Logic for within-group filtering
  keep_by_group <- rep(FALSE, nrow(TPM_matrix))
  for (grp in treatment_groups) {
    idx <- which(dge_obj$samples$treatment == grp)
    keep_grp <- rowSums(TPM_matrix[, idx, drop = FALSE] > TPM_cutoff) >= min_samples
    keep_by_group <- keep_by_group | keep_grp
  }
  # Logic for total expression across all samples
  keep_by_total <- rowSums(TPM_matrix) > total_cutoff
  # Combine: Must pass the group filter AND the total expression filter
  keep <- keep_by_group & keep_by_total
  return(keep)
}

# Subset samples first
dge_gene_gonad <- dge_gene[, dge_gene$samples$tissue == "Gonad"]
dge_gene_tail  <- dge_gene[, dge_gene$samples$tissue == "Tail"]    

# Drop unused factor levels
dge_gene_gonad$samples <- droplevels(dge_gene_gonad$samples)
dge_gene_tail$samples <- droplevels(dge_gene_tail$samples)         
# Apply treatment-aware filtering
keep_g <- filter_by_treatment(dge_gene_gonad)
keep_t <- filter_by_treatment(dge_gene_tail)                         

table(keep_g)
table(keep_t)

# Apply filters
dge_gene_gonad <- dge_gene_gonad[keep_g, , keep.lib.sizes = FALSE]
dge_gene_tail  <- dge_gene_tail[keep_t, , keep.lib.sizes = FALSE]  

# don't save for now, we will further filter these objects below

# Filter transcript-level count matrix

# Match gene-level keep vector to transcript-level GeneID
keep_tran_4 <- keep[dge$genes$GeneID]

# Filter dge object so that it only contains genes passing expression filter
dge_tran_4 <- dge[keep_tran_4, , keep.lib.sizes = FALSE]

# Remove transcripts with zero counts across all samples
nonzero_rows <- rowSums(dge_tran_4$counts) > 0
dge_tran_4 <- dge_tran_4[nonzero_rows, , keep.lib.sizes = FALSE]
dim(dge_tran_4)

# Normalize
dge_tran_4_norm <- calcNormFactors(dge_tran_4)

# Export count matrix and TMM-CPM matrix
tran_matrix_4 <- as.data.frame(dge_tran_4_norm$counts)
tran_matrix_CPM4 <- as.data.frame(cpm(dge_tran_4_norm, normalized.lib.sizes = TRUE), log = F)

# Add tran column
tran_matrix_4$transcript <- rownames(tran_matrix_4)
tran_matrix_CPM4$transcript <- rownames(tran_matrix_CPM4)

# Add gene column
tran_matrix_4$gene <- t2g$gene[match(tran_matrix_4$transcript, t2g$transcript)]
tran_matrix_CPM4$gene <- t2g$gene[match(tran_matrix_CPM4$transcript, t2g$transcript)]

# # Save
# write.csv(tran_matrix_4, 
#           here::here("data",
#                      "isoseq",
#                      "10_edgeR_filterexp",
#                      "01_all_tran_countmatrix_short.csv"))
# write.csv(tran_matrix_CPM4,
#           here::here("data",
#                      "isoseq",
#                      "10_edgeR_filterexp",
#                      "01_all_tran_TMMTPMmatrix_short.csv"))

# Look at isoform CPM distribution

# Calculate log2 CPM with prior.count = 1 (do this before normalizing)
logCPM <- cpm(dge, log=TRUE, prior.count=0.1)

# Convert to data frame and add gene IDs
df <- as.data.frame(logCPM)
df$TranscriptID <- rownames(logCPM)

# Gather into long format: columns = GeneID, Sample, logCPM
df_short <- df %>%
  pivot_longer(cols = -TranscriptID, names_to = "Sample", values_to = "logCPM")

# Plot
ggplot(df_short, aes(x = logCPM, color = Sample)) +
  geom_density() +
  labs(title = "Distribution of log2 CPM values per sample",
       x = "Log2 CPM",
       y = "Density") +
  theme_minimal()

# Save 
ggsave(here::here("figures", 
                  "05_SRP_histogram_isoformcpm_short.png"), 
       width = 12, 
       height = 10, 
       dpi = 300)

# Filter tissue-specific transcript matrices

# Split transcript-level into tissue specific objects as well
dge$samples$group <- factor(1)
dge_tran_gonad <- dge[, dge$samples$tissue == "Gonad"]
dge_tran_gonad$samples <- droplevels(dge_tran_gonad$samples)
dge_tran_tail<- dge[, dge$samples$tissue == "Tail"]
dge_tran_tail$samples <- droplevels(dge_tran_tail$samples)

# Match gene-level keep vector to transcript-level GeneID
keep_tran_g <- keep_g[dge_tran_gonad$genes$GeneID]
keep_tran_t <- keep_t[dge_tran_tail$genes$GeneID]

# Filter dge object so that it only contains genes passing expression filter
dge_tran_gonad <- dge_tran_gonad[keep_tran_g, , keep.lib.sizes = FALSE]
dge_tran_tail <- dge_tran_tail[keep_tran_t, , keep.lib.sizes = FALSE]

# Remove transcripts with zero counts across all samples
nonzero_rows <- rowSums(dge_tran_gonad$counts) > 0
dge_tran_gonad <- dge_tran_gonad[nonzero_rows, , keep.lib.sizes = FALSE]
nonzero_rows <- rowSums(dge_tran_tail$counts) > 0
dge_tran_tail <- dge_tran_tail[nonzero_rows, , keep.lib.sizes = FALSE]

# Check
dim(dge_tran_gonad)
# 54463    12
dim(dge_tran_tail)
# 40036    12

# Normalize using the TMM method in edgeR, then filter based on own criteria
# "TMM normalizes the library sizes to produce effective library sizes. 
# CPM values are counts normalized by the effective library sizes."
dge_tran_gonad <- calcNormFactors(dge_tran_gonad)
dge_tran_tail <- calcNormFactors(dge_tran_tail)

# Export count matrix and TMM-CPM matrix
tran_matrix_gonad <- as.data.frame(dge_tran_gonad$counts)
tran_matrix_CPMgonad <- as.data.frame(cpm(dge_tran_gonad, normalized.lib.sizes = TRUE), log = F)
tran_matrix_tail <- as.data.frame(dge_tran_tail$counts)
tran_matrix_CPMtail <- as.data.frame(cpm(dge_tran_tail, normalized.lib.sizes = TRUE), log = F)

# Add tran column
tran_matrix_gonad$transcript <- rownames(tran_matrix_gonad)
tran_matrix_CPMgonad$transcript <- rownames(tran_matrix_CPMgonad)
tran_matrix_tail$transcript <- rownames(tran_matrix_tail)
tran_matrix_CPMtail$transcript <- rownames(tran_matrix_CPMtail)

# Add gene column
tran_matrix_gonad$gene    <- t2g$gene[match(tran_matrix_gonad$transcript, t2g$transcript)]
tran_matrix_CPMgonad$gene <- t2g$gene[match(tran_matrix_CPMgonad$transcript, t2g$transcript)]
tran_matrix_tail$gene     <- t2g$gene[match(tran_matrix_tail$transcript, t2g$transcript)]
tran_matrix_CPMtail$gene  <- t2g$gene[match(tran_matrix_CPMtail$transcript, t2g$transcript)]

# Plot MDS (this will automatically do a log cpm transform)
# plotMDS(dge_gene_gonad,
#         col = as.numeric(dge_gene_gonad$samples$sex),        # color by sex
#         labels = dge_gene_gonad$samples$sample)              # label with sample name
# plotMDS(dge_gene_tail,
#         col = as.numeric(dge_gene_tail$samples$sex),        # color by sex
#         labels = dge_gene_tail$samples$sample)              # label with sample name
# plotMDS(dge_tran_gonad,
#         col = as.numeric(dge_tran_gonad$samples$sex),        # color by sex
#         labels = dge_tran_gonad$samples$sample)              # label with sample name
# plotMDS(dge_tran_tail,
#         col = as.numeric(dge_tran_tail$samples$sex),        # color by sex
#         labels = dge_tran_tail$samples$sample)              # label with sample name

# Save
write.csv(tran_matrix_gonad, 
          here::here("data",
                     "02_build_shortread_transcriptome",
                     "05_edgeR",
                     "01_gonadtran_countmatrix_short.csv"))
write.csv(tran_matrix_CPMgonad,
          here::here("data",
                     "02_build_shortread_transcriptome",
                     "05_edgeR",
                     "01_gonadtran_TMMTPMmatrix_short.csv"))
write.csv(tran_matrix_tail, 
          here::here("data",
                     "02_build_shortread_transcriptome",
                     "05_edgeR",
                     "02_tailtran_countmatrix_short.csv"))
write.csv(tran_matrix_CPMtail,
          here::here("data",
                     "02_build_shortread_transcriptome",
                     "05_edgeR",
                     "02_tailtran_TMMTPMmatrix_short.csv"))

# Wrangle treatment level count matrix

# Pivot_shorter so that you have a sample column
tran_matrix_CPM4_short <- tran_matrix_CPM4 %>%
  pivot_longer(cols = matches("^R[FM][GT][0-9]$"),  # sample columns (RFG1 etc.)
               names_to = "sample",         # Create 'sample' column
               values_to = "isoform_cpm")   # Create 'cpm' column with values

# Aggregate transcript counts by gene
collapsed_gene_matrix <- tran_matrix_CPM4_short %>%
  group_by(gene, sample) %>%
  summarise(gene_cpm = sum(isoform_cpm, na.rm = TRUE))

# Merge with metadata
collapsed_gene_matrix <- metadata %>%
  dplyr::select(sample_id, sex, tissue) %>%
  right_join(collapsed_gene_matrix, by = c("sample_id" = "sample"))

# Join with isoform data
combined_data <- collapsed_gene_matrix %>%
  left_join(tran_matrix_CPM4_short,
            by = c("sample_id" = "sample", 
                   "gene" = "gene"))

################################################################################
### (2) Explore expression
################################################################################
# 
# # Plot # isoform histogram-
# 
# # Remove any isoforms with no expression
# 
# # Find # isoforms per gene
# cpm_v_numisoforms <- combined_data %>%
#   group_by(gene) %>%
#   summarise(
#     num_isoforms = n_distinct(transcript) 
#   ) %>%
#   distinct()
# 
# # Histogram isoforms per gene
# cpm_v_numisoforms %>%
#   ggplot(aes(x = num_isoforms)) +
#   geom_histogram(binwidth = 1, boundary = 0.5, fill = "#a45851", color = "black") +
#   scale_x_continuous(breaks = scales::pretty_breaks()) +
#   labs(
#     x = "Number of Isoforms",
#     y = "Count of Genes"
#   ) +
#   theme_minimal()
# 
# # Save
# ggsave(
#   filename = here::here("figures",
#                         "05_SRP_histo_numisoforms_short.png"),
#   width = 8,           # width in inches
#   height = 6,          # height in inches
#   dpi = 300,           # resolution (use 300 for publication-quality)
# )
# 
# # Sum gene_cpm across all samples
# combined_data <- combined_data %>%
#   dplyr::select(sample_id, gene, gene_cpm) %>%  # Select only relevant columns
#   distinct() %>%  # Keep unique combinations of sample and gene
#   group_by(gene) %>%  # Group by gene to sum across all samples
#   summarise(total_gene_cpm = sum(gene_cpm, na.rm = TRUE)) %>%  # Sum gene_cpm across samples for each gene
#   ungroup() %>% # Remove grouping
#   right_join(combined_data, 
#              by = c("gene"))
# 
# # Sum isoform_cpm across all samples
# combined_data <- combined_data %>%
#   group_by(transcript) %>%
#   mutate(total_isoform_cpm = sum(isoform_cpm, na.rm = TRUE)) %>%
#   ungroup()
# 
# # Calculate relative abundance per isoform
# combined_data <- combined_data %>%
#   mutate(total_relabund_isoform = total_isoform_cpm / total_gene_cpm,
#          sample_relabund_isoform = isoform_cpm / gene_cpm)
# 
# # Plot relative abundance of all isoforms
# 
# # Summarize the data to ensure we don't double-count isoforms
# combined_data_unique <- combined_data %>%
#   dplyr::select(gene, transcript, total_gene_cpm, total_isoform_cpm, total_relabund_isoform) %>%
#   distinct()
# 
# # Plot
# ggplot(combined_data_unique, aes(x = total_relabund_isoform)) +
#   geom_histogram(breaks = seq(0, 1, by = .025), fill = "#A8B28D", color = "black") +
#   labs(
#     x = "Isoform Relative Abundance",
#     y = "Count"
#   ) +
#   theme_minimal()
# 
# # Save
# ggsave(
#   filename = here::here("figures",
#                         "05_isoform_relabund_counts_short.png"),
#   width = 8,           # width in inches
#   height = 6,          # height in inches
#   dpi = 300,           # resolution (use 300 for publication-quality)
# )
# 
# # Plot using percentage to match Benetiere et al. 2024
# ggplot(combined_data_unique, aes(x = total_relabund_isoform)) +
#   geom_histogram(aes(y = after_stat(count / sum(count))), 
#                  breaks = seq(0, 1, by = .05), fill = "#A8B28D", color = "black", boundary = 0) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   labs(
#     x = "Isoform relative abundance",
#     y = "Percent"
#   ) +
#   theme_minimal()
# 
# # Save
# ggsave(
#   filename = here::here("figures",
#                         "05_isoform_relabund_percent_short.png"),
#   width = 8,           # width in inches
#   height = 6,          # height in inches
#   dpi = 300,           # resolution (use 300 for publication-quality)
# )
# 
# # Calculate quartile counts again
# combined_data_unique %>%
#   mutate(ratio_quartile = cut(
#     total_relabund_isoform,
#     breaks = c(0, 0.25, 0.5, 0.75, 1),
#     include.lowest = TRUE,
#     labels = c("0–0.25", "0.25–0.5", "0.5–0.75", "0.75–1")
#   )) %>%
#   dplyr::count(ratio_quartile, name = "num_isoforms")
# 
# # # Save
# # write.csv(combined_data_unique, 
# #           here::here("data",
# #                      "isoseq",
# #                      "10_edgeR_filterexp",
# #                      "01_all_collapsed_relabund_short.csv"))

# Save edgeR objects that have been filtered

### After determining which isoforms are effective in script 5, 
### we will further filter these and then do DGE, DTE, DTU analyses

# Sanity check
dim(dge_gene_gonad)
dim(dge_gene_tail)
dim(dge_tran_gonad)
dim(dge_tran_tail)

# Save edgeR objects to a .RData file
save(dge_gene_gonad, file = here::here("data",
                                       "02_build_shortread_transcriptome",
                                       "05_edgeR",
                                       "01_eR_genegonad_short.RData"))
save(dge_gene_tail, file = here::here("data",
                                      "02_build_shortread_transcriptome",
                                      "05_edgeR",
                                      "02_eR_genetail_short.RData"))
save(dge_tran_gonad, file = here::here("data",
                                       "02_build_shortread_transcriptome",
                                       "05_edgeR",
                                       "01_eR_trangonad_short.RData"))
save(dge_tran_tail, file = here::here("data",
                                      "02_build_shortread_transcriptome",
                                      "05_edgeR",
                                      "02_eR_trantail_short.RData"))

################################################################################
### (3) Filter over files by filter outputs
################################################################################

# Get the list of transcripts that passed in tail and gonad to keep
keep_transcripts_union <- unique(c(
  tran_matrix_CPMgonad$transcript, 
  tran_matrix_CPMtail$transcript
))

# Filter GRanges object 
gtf_filtered <- gtf[gtf$transcript_id %in% keep_transcripts_union]

# Check counts
cat("Genes:", length(unique(gtf_filtered$gene_id)), "\n")
cat("Transcripts:", length(unique(gtf_filtered$transcript_id)), "\n")

# Export as valid GTF (will automatically use the proper 9-column format)
rtracklayer::export(
  gtf_filtered,
  here::here("data", 
       "GTF", 
       "06_shortreadpipe_filtered_by_exp_short.gtf"),
  format = "gtf"
)
