## Step 16: edgeR expression filter on the long-read kallisto counts
## Last modified by Linley Sherin 2026-10-27
#
# FILTER NOTICE: Filter genes without >1TPM in 3/6 samples of a sex
#
# Runs locally.
# (1) Load raw count matrices (kallisto) - includes overdispersion correction
# (2) Create DGE objects, inspect CPM distribution and determine CPM cut-off
# (3) Subset DGE objects by tissue and apply filter at gene level
# (4) Save filtered DGE objects for DE
# (5) Filter over files by filter outputs
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 01_build_longread_transcriptome/16_edgeR_filterexp.R
#
# Reads:  data/01_build_longread_transcriptome/15_kallisto/kallisto_long/flnc_<sample>.fastq/quant_unfiltered/abundance_1/
#         data/01_build_longread_transcriptome/metadata.csv
#         data/GTF/03_nonredundant_transcriptome.gtf
#         data/01_build_longread_transcriptome/12_sqanti_filter/SQANTI_filtered_RulesFilter_result_classification.txt
#         data/01_build_longread_transcriptome/14_remove_redundancy/rosetta_stone.csv
# Writes: data/01_build_longread_transcriptome/16_edgeR_filterexp/
#           dispersion_output.csv
#           01_gonadtran_*.csv, 02_tailtran_*.csv     count and TMM-TPM matrices
#           01_eR_gene/trangonad.RData, 02_eR_gene/trantail.RData   filtered DGE objects
#           SQANTI_classification_filtered.csv
#           04_filtered_by_exp_long.gtf
#         data/GTF/04_filtered_by_exp_long.gtf
#         figures/16_histogram_geneTPM_long.png, 16_histogram_isoformTPM_long.png
#
# Overwrites its outputs every time it runs.

library(tidyr)
library(ggplot2)
library(dplyr)
library(edgeR)
library(rtracklayer)

################################################################################
################################################################################
### LONG READS
################################################################################
################################################################################

################################################################################
### (1) Load raw count matrices (kallisto)
################################################################################

# edgeR wants the directory where the abundance.h5 files are kept as an input

# List parent directory
base_dir <- here::here("data",
                       "01_build_longread_transcriptome",
                       "15_kallisto",
                       "kallisto_long")

# Get all sample directory names
sample_dirs <- list.dirs(base_dir, 
                         recursive = FALSE)

# Build the full paths to abundance_1 folders
abundance_dirs <- file.path(sample_dirs, 
                            "quant_unfiltered", 
                            "abundance_1")

# Also load metadata for batch effects
# Contains: cell, sample_id, sex, tissue, selection, read depth (from flnc.bam)
metadata <- read.csv(here::here("data",
                                "01_build_longread_transcriptome",
                                "metadata.csv"))

# Read in the data using catchKallisto
y <- catchKallisto(abundance_dirs)

# Load gtf
gtf <- import(here::here("data",
                         "GTF",
                         "03_nonredundant_transcriptome.gtf"))
gtf_df <- as.data.frame(gtf)

### FILTER NOTICE: Remove all non-autosomal chromosomes
# Remove non-autosomal genes (Chr 21 and unplaced scaffolds)
gtf_df_filtered <- gtf_df[!grepl("^NW", as.character(gtf_df$seqnames)), ] %>%
  filter(seqnames != "NC_024342.1") %>% # remove transcripts on 12
  filter(type == "transcript") # only keep transcript rows

# Summary
cat("Number of genes before filtering:", length(unique(gtf_df_filtered$gene_id)), "\n")
cat("Number of genes before filtering:", length(unique(gtf_df_filtered$transcript_id)), "\n")

# Filter
y[c("counts", "annotation")] <- lapply(
  y[c("counts", "annotation")],
  function(x) x[rownames(x) %in% unique(gtf_df_filtered$transcript_id), , drop = FALSE]
)

# Create a data frame with gene names and dispersions
dispersion_df <- data.frame(
  transcript = rownames(y$annotation),     # transcript ids from counts matrix
  dispersion = y$annotation$Overdispersion    # tagwise dispersions (transcript level)
)

# Save
write.csv(dispersion_df,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "dispersion_output.csv"))

### Create metadata

# Your vector of all sample names
all_samples <- c("RFG1", "RFG2", "RFG3", "RFG4", "RFG5", "RFG6",
                 "RFT1", "RFT2", "RFT3", "RFT4", "RFT5", "RFT6",
                 "RMG1", "RMG2", "RMG3", "RMG4", "RMG5", "RMG6",
                 "RMT1", "RMT2", "RMT3", "RMT4", "RMT5", "RMT6")

# Combine into a metadata data frame
sample_info <- data.frame(
  sample = all_samples,
  sex = sex <- factor(ifelse(grepl("F", all_samples), "Female", "Male")),
  tissue = tissue <- factor(ifelse(grepl("T", all_samples), "Tail", "Gonad")),
  treatment = factor(paste0(sex, "_", tissue))
)

# Add batch effects from metadata (cell) and selection column
sample_info <- sample_info %>%
  left_join(metadata %>% 
              dplyr::select(sample_id, selection, cell),
            by = c("sample" = "sample_id"))

################################################################################
### (2) Create DGE objects, inspect CPM distribution and determine CPM cut-off
################################################################################

### Create transcript-level DGE objects

# Use dividied count strategy to control for Overdispersion at transcript level
dge <- DGEList(counts = y$counts / y$annotation$Overdispersion,
               genes = y$annotation,
               samples = sample_info,
               group = sample_info$treatment)

### Add sample info to dge object (fixing directory names)
rownames(sample_info) <- sample_info$sample
# Replace column names with the part after flnc_ and before the next "\\"
colnames(dge) <- sub(".*flnc_(R[FM][GT][0-9]).*", "\\1", colnames(dge))
# Check the new column names
print(colnames(dge))
all(rownames(sample_info) == colnames(dge))  # should return TRUE
dge$samples <- cbind(sample_info[colnames(dge), ], 
                     dge$samples[, c("lib.size", "norm.factors")])

### Collapse to gene-level

# Add gene id and transcript id to dge
dge$genes$TranscriptID <- rownames(dge)
dge$genes$GeneID <- sub("\\.[^.]+$", "", rownames(dge))

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

# Remove columns we no longer need
dge_gene$genes$TranscriptID <- NULL
dge_gene$samples$norm.factors.1 <- NULL
dge_gene$samples$lib.size.1 <- NULL

### Look at gene-level TPM distribution

counts <- data.frame(dge_gene)
sample_totals <- colSums(counts[, grep("^R(FG|FT|MG|MT)[0-9]+$", names(counts))])
data.frame(sample_totals)

# Calculate log2 TPM with prior.count = 1 (do this before normalizing)
logTPM <- cpm(dge_gene, log=TRUE, prior.count=0.1)

# Convert to data frame and add gene IDs
df <- as.data.frame(logTPM)
df$GeneID <- rownames(logTPM)

# Gather into long format: columns = GeneID, Sample, logTPM
df_long <- df %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "logTPM")

# Plot
ggplot(df_long, aes(x = logTPM, color = Sample)) +
  geom_density() +
  labs(title = "Distribution of log2 TPM values per sample",
       x = "Log2 TPM",
       y = "Density") +
  theme_minimal()

# Save 
ggsave(here::here("figures", 
                  "16_histogram_geneTPM_long.png"), 
       width = 12, 
       height = 10, 
       dpi = 300)

### USE 1 AS A TPM CUT-OFF

################################################################################
### (3) Subset DGE objects by tissue and apply filter at gene level
################################################################################

### Subset DGE into each pairwise comparison in a new objects and filter

### Originally I thought I needed to filter at the dataset level (collapsed tissues)
### and also tissue level, but now I am just going to do tissue level
### because I think the other way complicates things too much

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
# dim(dge_gene_gonad)
# [1] 16630    12

# Drop unused factor levels
dge_gene_gonad$samples <- droplevels(dge_gene_gonad$samples)
dge_gene_tail$samples  <- droplevels(dge_gene_tail$samples)

# Apply treatment-aware filtering
keep_g <- filter_by_treatment(dge_gene_gonad)
keep_t <- filter_by_treatment(dge_gene_tail)
table(keep_g)
table(keep_t)

# Apply filters
dge_gene_gonad <- dge_gene_gonad[keep_g, , keep.lib.sizes = FALSE]
dge_gene_tail  <- dge_gene_tail[keep_t, , keep.lib.sizes = FALSE]

# Normalize using the TMM method in edgeR, then filter based on own criteria
dge_gene_gonad <- calcNormFactors(dge_gene_gonad)
dge_gene_tail <- calcNormFactors(dge_gene_tail)

# Check
dge_gene_gonad$samples$norm.factors
dge_gene_tail$samples$norm.factors

# Export count matrix and TMM-TPM matrix
gene_matrix_gonad <- as.data.frame(dge_gene_gonad$counts)
gene_matrix_TPMgonad <- as.data.frame(cpm(dge_gene_gonad, normalized.lib.sizes = TRUE), log = F)
gene_matrix_tail <- as.data.frame(dge_gene_tail$counts)
gene_matrix_TPMtail <- as.data.frame(cpm(dge_gene_tail, normalized.lib.sizes = TRUE), log = F)

### Look at isoform TPM distribution

# Calculate log2 TPM with prior.count = 1 (do this before normalizing)
logTPM <- cpm(dge, log=TRUE, prior.count=0.1)

# Convert to data frame and add gene IDs
df <- as.data.frame(logTPM)
df$TranscriptID <- rownames(logTPM)

# Gather into long format: columns = GeneID, Sample, logTPM
df_long <- df %>%
  pivot_longer(cols = -TranscriptID, names_to = "Sample", values_to = "logTPM")

# Plot
ggplot(df_long, aes(x = logTPM, color = Sample)) +
  geom_density() +
  labs(title = "Distribution of log2 TPM values per sample",
       x = "Log2 TPM",
       y = "Density") +
  theme_minimal()

# Save 
ggsave(here::here("figures", 
                  "16_histogram_isoformTPM_long.png"), 
       width = 12, 
       height = 10, 
       dpi = 300)

### Filter tissue-specific transcript matrices

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

# Normalize using the TMM method in edgeR, then filter based on own criteria
dge_tran_gonad <- calcNormFactors(dge_tran_gonad)
dge_tran_tail <- calcNormFactors(dge_tran_tail)

# Check
dge_tran_gonad$samples$norm.factors
dge_tran_tail$samples$norm.factors

# Export count matrix and TMM-TPM matrix
tran_matrix_gonad <- as.data.frame(dge_tran_gonad$counts)
tran_matrix_TPMgonad <- as.data.frame(cpm(dge_tran_gonad, normalized.lib.sizes = TRUE), log = F)
tran_matrix_tail <- as.data.frame(dge_tran_tail$counts)
tran_matrix_TPMtail <- as.data.frame(cpm(dge_tran_tail, normalized.lib.sizes = TRUE), log = F)

# Add tran column
tran_matrix_gonad$transcript <- rownames(tran_matrix_gonad)
tran_matrix_TPMgonad$transcript <- rownames(tran_matrix_TPMgonad)
tran_matrix_tail$transcript <- rownames(tran_matrix_tail)
tran_matrix_TPMtail$transcript <- rownames(tran_matrix_TPMtail)

# Add gene column
tran_matrix_gonad$gene <- sub("\\.[^.]+$", "", tran_matrix_gonad$transcript)
tran_matrix_TPMgonad$gene <- sub("\\.[^.]+$", "", tran_matrix_TPMgonad$transcript)
tran_matrix_tail$gene <- sub("\\.[^.]+$", "", tran_matrix_tail$transcript)
tran_matrix_TPMtail$gene <- sub("\\.[^.]+$", "", tran_matrix_TPMtail$transcript)

# Save
write.csv(tran_matrix_gonad, 
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "01_gonadtran_countmatrix.csv"))
write.csv(tran_matrix_TPMgonad,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "01_gonadtran_TMMTPMmatrix.csv"),
          row.names = F)
write.csv(tran_matrix_tail, 
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "02_tailtran_countmatrix.csv"))
write.csv(tran_matrix_TPMtail,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "02_tailtran_TMMTPMmatrix.csv"))

# Plot MDS (this will automatically do a log TPM transform)
plotMDS(dge_gene_gonad,
         col = as.numeric(dge_gene_gonad$samples$sex),        # color by sex
         labels = dge_gene_gonad$samples$sample)              # label with sample name
plotMDS(dge_gene_tail,
         col = as.numeric(dge_gene_tail$samples$sex),        # color by sex
         labels = dge_gene_tail$samples$sample)              # label with sample name
plotMDS(dge_tran_gonad,
         col = as.numeric(dge_tran_gonad$samples$sex),        # color by sex
         labels = dge_tran_gonad$samples$sample)              # label with sample name
plotMDS(dge_tran_tail,
         col = as.numeric(dge_tran_tail$samples$sex),        # color by sex
         labels = dge_tran_tail$samples$sample)              # label with sample name

################################################################################
### (4) Save filtered DGE objects for DE
################################################################################

### Save edgeR objects that have been filtered

### After determining which isoforms are effective in script 5, 
### we will further filter these and then do DGE, DTE, DTU analyses

# Sanity check
dim(dge_gene_gonad)
dim(dge_gene_tail)
dim(dge_tran_gonad)
dim(dge_tran_tail)

# Save edgeR objects to a _long.RData file
save(dge_gene_gonad, file = here::here("data",
                                       "01_build_longread_transcriptome",
                                       "16_edgeR_filterexp",
                                       "01_eR_genegonad.RData"))
save(dge_gene_tail, file = here::here("data",
                                       "01_build_longread_transcriptome",
                                       "16_edgeR_filterexp",
                                       "02_eR_genetail.RData"))
save(dge_tran_gonad, file = here::here("data",
                                       "01_build_longread_transcriptome",
                                       "16_edgeR_filterexp",
                                       "01_eR_trangonad.RData"))
save(dge_tran_tail, file = here::here("data",
                                       "01_build_longread_transcriptome",
                                       "16_edgeR_filterexp",
                                       "02_eR_trantail.RData"))

################################################################################
### (5) Filter over files by filter outputs
################################################################################

# Get the list of transcripts that passed in tail and gonad to keep
keep_transcripts_union <- unique(c(
  tran_matrix_TPMgonad$transcript, 
  tran_matrix_TPMtail$transcript
))

# Classification.txt output from SQANTI filter
class <- read.delim2(here::here("data",
                                "01_build_longread_transcriptome",
                                "12_sqanti_filter",
                                "SQANTI_filtered_RulesFilter_result_classification.txt"),
                  header = T)
class <- class[class$filter_result == "Isoform", ] # Filter out isoforms that failed SQANTI filter

# Read in gene reassignment map to adjust files
map <- read.csv(here::here("data",
                           "01_build_longread_transcriptome",
                           "14_remove_redundancy",
                           "rosetta_stone.csv"))

# Merge on transcript_id (same as class$isoform)
class <- left_join(class, map, by = c("isoform" = "transcript_id"))

# Rename new columns
class <- class %>%
  dplyr::rename(
    gene = new_gene_id,
    transcript = new_transcript_id
  )

# Remove columns where all values are NA
class <- class[, colSums(!is.na(class)) > 0]

# Remove FL.* columns using dplyr::select and tidyselect
class <- class %>% dplyr::select(-starts_with("FL."))

# Now filter to match the dfs after expression filter (match keep_tran_4)
# Keep only rows where transcript is in the DGEList
class <- class %>% 
  filter(transcript %in% keep_transcripts_union)

# Save class
write.csv(class,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "SQANTI_classification_filtered.csv"),
          row.names = F)

# # GTF output from redundancy script
# gtf <- rtracklayer::import(here::here("data",
#                                       "GTF",
#                                       "03_nonredundant_transcriptome.gtf"))

# Filter GRanges object 
gtf_filtered <- gtf[gtf$transcript_id %in% keep_transcripts_union]

# Export as valid GTF (will automatically use the proper 9-column format)
rtracklayer::export(
  gtf_filtered,
  here::here("data", 
       "GTF", 
       "04_filtered_by_exp_long.gtf"),
  format = "gtf"
)

# Keep a copy in this step's folder too
file.copy(here::here("data", "GTF", "04_filtered_by_exp_long.gtf"),
          here::here("data", "01_build_longread_transcriptome", "16_edgeR_filterexp", "04_filtered_by_exp_long.gtf"),
          overwrite = TRUE)

# Summary
cat("Number of genes after filtering:", length(unique(gtf_filtered$gene_id)), "\n")
cat("Number of transcripts after filtering:", length(unique(gtf_filtered$transcript_id)), "\n")
