# Subsample count matrixes to see whether we reached isoform saturation

# Saturation curves

library(data.table)
library(ggplot2)
library(tidyverse)
library(here)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(stringr)
library(scales)
library(colorspace)

# Load data -------------------------------------------------------------------

# Load TPM matrices for gonad and tail
count_gonad <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "16_edgeR_filterexp",
                                 "01_gonadtran_countmatrix.csv"))[,-1] %>%
  rename_with(.fn = ~ sub("^X", "", .), .cols = 1:14)

count_tail <- read.csv(here::here("data",
                                "01_build_longread_transcriptome",
                                "16_edgeR_filterexp",
                                "02_tailtran_countmatrix.csv"))[,-1] %>%
  rename_with(.fn = ~ sub("^X", "", .), .cols = 1:14)

# Combine them
count_matrix <- full_join(count_gonad, count_tail)

# Convert NAs to 0s
count_matrix[is.na(count_matrix)] <- 0

# Identify the columns 'gene' and 'transcript'
gene_col <- "gene"
transcript_col <- "transcript"

# Get the current column order
current_columns <- colnames(count_matrix)

# Remove 'gene' and 'transcript' from the list of columns
remaining_columns <- setdiff(current_columns, c(gene_col, transcript_col))

# Reorder the columns so 'gene' and 'transcript' come at the end
new_column_order <- c(remaining_columns, gene_col, transcript_col)

# Set the new column order
setcolorder(count_matrix, new_column_order)

# Define colours ----------------------------------------------------------

# Define your colors for each sex_tissue group
custom_colors <- c(
  "F_G" = "#C08081",  
  "F_T" = "#D9B44A",  
  "M_G" = "#7E9E7E",  
  "M_T" = "#8C7A6B"   
)

# Calculate the mapping rate correction for each sample ------------------------

# Contains: cell, sample_id, sex, tissue, selection, read depth (from flnc.bam)
metadata <- read.csv(here::here("data",
                                "isoseq",
                                "metadata.csv"))

# Calculate mapped_readdepth
mapped_readdepth <- count_matrix %>%
  dplyr::select(-transcript, -gene) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "sample_id", values_to = "mapped_readdepth") %>%
  mutate(sample_id = str_remove(sample_id, "_isoform_counts"))  # To match metadata

metadata$flnc_readdepth <- as.numeric(metadata$flnc_readdepth)

# Join mapped_readdepth with metadata to calculate mapping rate
metadata <- metadata %>%
  left_join(mapped_readdepth, by = "sample_id") %>%
  dplyr::mutate(mapping_rate = mapped_readdepth / flnc_readdepth)

# Calculate mapping rate correction across all samples ###########--------------

# We want to correct the read depth #, which currently represents the reads 
# mapped by minimap2 and quantified by bambu so that it represents actual
# sample effort (i.e. the # reads we payed for and sequenced)

# We have to account for two things
# (1) reads lost from raw -> end of pacbio workflow (7.25%)
# (2) reads lost from pacbio workflow -> after minimap2
raw_reads <- 195014391 * 1.0725
flnc_reads <- 195014391 # this is the input for minimap2

setDT(count_matrix)  # Convert to data.table

# Reads in count_matrix represents reads after kallisto and filtering by low TPM
count_matrix_prune <- count_matrix[, !(names(count_matrix) %in% c("transcript", "gene")), with = FALSE]

# Convert to numeric matrix
# Identify columns whose names end with 'p'
p_columns <- grep("P$", colnames(count_matrix), value = TRUE)

# Convert these columns to numeric
count_matrix[, (p_columns) := lapply(.SD, as.numeric), .SDcols = p_columns]

# Sum all reads
after_reads <- sum(as.matrix(count_matrix_prune), na.rm = TRUE)

# So now we can calculate the percent of reads that made it from raw to the final step (count matrix)
mapping_rate <- after_reads / raw_reads

# Subsample count matrix by sample ---------------------------------------------

# Define the sampling depths (this will sample every 10% of reads)
fractions <- seq(0.1, 1.0, by = 0.05)

# Pre-allocate an empty data.table to fill by rows
saturation_dt <- data.table(
  sample = character(),
  fraction = numeric(),
  total_reads = integer(),
  detected_transcripts = integer(),
  detected_genes = integer()
)

# transcript IDs and Gene IDs
transcript_ids <- count_matrix$transcript
gene_ids <- count_matrix$gene

# Loop over samples and fractions efficiently
for (sample in colnames(count_matrix)[1:24]) {
  
  sample_counts <- count_matrix[[sample]]
  sample_counts[is.na(sample_counts)] <- 0
  
  total_reads <- sum(sample_counts)
  
  # Skip the iteration if total_reads == 0
  if (total_reads == 0) {
    next
  }
  
  probs <- sample_counts / total_reads
  
  temp_res <- data.table(
    sample = sample,
    fraction = fractions,
    total_reads = round(total_reads * fractions),
    detected_transcripts = integer(length(fractions)),
    detected_genes = integer(length(fractions))
  )
  
  temp_res[, c("detected_transcripts", "detected_genes") := {
    results <- lapply(total_reads, function(n_reads) {
      subsampled <- rmultinom(1, size = n_reads, prob = probs)
      detected_transcripts <- sum(subsampled > 0)
      detected_genes <- length(unique(gene_ids[subsampled > 0]))
      list(detected_transcripts, detected_genes)
    })
    list(sapply(results, `[[`, 1), sapply(results, `[[`, 2))
  }]
  
  saturation_dt <- rbind(saturation_dt, temp_res)
}

# Change sample id
saturation_dt[, sample := stringr::str_remove(sample, "_transcript_counts")]

# Merge saturation results with sample-specific mapping rates
saturation_dt <- merge(
  saturation_dt,
  metadata[, c("sample_id", "mapping_rate")],  # Base R style
  by.x = "sample", 
  by.y = "sample_id",
  all.x = TRUE
)

# Create corrected read depth per sample using sample-specific mapping rate
saturation_dt[, total_reads_corrected := round(total_reads / mapping_rate)]

# Define groups and plot at the sample level -----------------------------------

# Extract groupings/treatments using regex groups
saturation_dt <- saturation_dt %>%
  mutate(
    # Extract with full match groups: number, sex, tissue
    match = str_match(sample, "([0-9]+)([FM])([GT])"),
    
    sample_num = as.numeric(match[, 2]),  # first group
    sex        = match[, 3],              # second group
    tissue     = match[, 4],              # third group
    
    # Determine selection
    selection = if_else(sample_num %% 2 == 1, "up", "down"),
    
    # Combine to final group
    treatment = paste0(sex, "_", tissue)
  ) %>%
  dplyr::select(-match)  # optional: drop intermediate match column


# Plot
ggplot(saturation_dt, aes(x = total_reads_corrected, y = detected_transcripts, group = sample, color = treatment)) +
  geom_line(alpha = 0.9) +       # individual lines for each sample, with some transparency
  geom_point(size = 1.5, alpha = 0.9) + # optional points for clarity
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 50000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Read Depth",
       y = "Detected Isoforms",
       color = "Group") +
  theme_minimal() +
  theme(legend.position = "right")

# Save
ggsave(filename = here::here("figures",
                             "10_satcurve_sample.png"),
       width = 6, 
       height = 4, 
       dpi = 300)

# Plot
ggplot(saturation_dt, aes(x = total_reads_corrected, y = detected_genes, group = sample, color = treatment)) +
  geom_line(alpha = 0.9) +       # individual lines for each sample, with some transparency
  geom_point(size = 1.5, alpha = 0.9) + # optional points for clarity
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = scales::comma, limits = c(5000, 20000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Read Depth",
       y = "Detected Genes",
       color = "Group") +
  theme_minimal() +
  theme(legend.position = "right")

# Save
ggsave(filename = here::here("figures",
                             "10_satcurve_sample_gene.png"), 
       width = 6, 
       height = 4, 
       dpi = 300)


# Plot again by proportion transcripts detected -----------------------------------
# 
# # Calculate total transcripts at fraction == 1 for each sample
# total_transcripts_per_sample <- saturation_dt[fraction == 1, .(total_transcripts = detected_transcripts), by = sample]
# 
# # Join back and calculate proportion
# saturation_dt <- merge(saturation_dt, total_transcripts_per_sample, by = "sample")
# saturation_dt[, prop_transcript := detected_transcripts / total_transcripts]
# 
# # Plot
# ggplot(saturation_dt, aes(x = total_reads_corrected, 
#                           y = prop_transcript, 
#                           group = sample, 
#                           color = treatment)) +
#   geom_line(alpha = 0.9) +
#   geom_point(size = 1.5, alpha = 0.9) +
#   scale_color_manual(values = custom_colors) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.5, 1)) +
#   scale_x_continuous(labels = scales::comma) +
#   labs(
#     x = "Read Depth",
#     y = "Proportion of transcripts Detected",
#     color = "Group"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# # Save
# ggsave(filename = here::here("figures", "10_satcurve_prop_sample.png"), 
#        width = 6, 
#        height = 4, 
#        dpi = 300)
# 

# Subsample count matrix by tissue ---------------------------------------------

# Convert count_matrix to data.table
setDT(count_matrix)

# Extract the sample columns (excluding transcript and gene)
sample_cols <- setdiff(colnames(count_matrix), c("transcript", "gene"))

# Create mapping of sample -> treatment
sample_info <- data.frame(sample = sample_cols) %>%
  mutate(
    # Extracting sample_num, sex, tissue, and selection
    match = str_match(sample, "^([0-9]+)([FM])([A-Z]+)"),
    sample_num = as.integer(match[, 2]),
    sex = if_else(match[, 3] == "F", "female", "male"),
    tissue = if_else(match[, 4] == "GP", "gonad", "tail"),
    selection = if_else(sample_num %% 2 == 1, "up", "down"),
    treatment = paste(sex, tissue, sep = "_")
  )

# Create collapsed count matrix by group
counts_collapsed <- data.table(transcript = count_matrix$transcript)

for (group in unique(sample_info$treatment)) {
  # Get all sample columns for this group
  group_samples <- sample_info$sample[sample_info$treatment == group]
  
  # Sum counts for the group
  counts_collapsed[[group]] <- rowSums(count_matrix[, ..group_samples], na.rm = TRUE)
}

# Create output table for saturation results
saturation_group <- data.table()

# Define your fractions ahead of time if not already defined
# e.g., fractions <- seq(0.1, 1.0, by = 0.1)

for (group in colnames(counts_collapsed)[-1]) {  # skip transcript
  group_counts <- counts_collapsed[[group]]
  group_counts[is.na(group_counts)] <- 0
  total_reads <- sum(group_counts)
  if (total_reads == 0) next
  
  probs <- group_counts / total_reads
  
  temp_res <- data.table(
    treatment = group,
    fraction = fractions,
    total_reads = round(total_reads * fractions),
    detected_transcripts = integer(length(fractions))
  )
  
  temp_res[, detected_transcripts := sapply(total_reads, function(n_reads) {
    subsampled <- rmultinom(1, size = n_reads, prob = probs)
    sum(subsampled > 0)
  })]
  
  saturation_group <- rbind(saturation_group, temp_res)
}

# Convert to data.table if not already
setDT(metadata)

# Define grouping label (same format as in saturation_group)
metadata[, treatment := paste(sex, tissue, sep = "_")]

# Compute mean mapping rate per group
group_mapping_rates <- metadata[, .(
  mapping_rate = mean(mapping_rate, na.rm = TRUE)
), by = treatment]

# Merge group mapping rates into saturation_group
saturation_group <- merge(
  saturation_group,
  group_mapping_rates,
  by = "treatment",
  all.x = TRUE
)

# Calculate corrected total reads
saturation_group[, total_reads_corrected := round(total_reads / mapping_rate)]

# Define groups and plot at the tissue level -----------------------------------

# Redefine to match the current set up
custom_colors <- c(
  "female_gonad" = "#C08081",  
  "female_tail" = "#D9B44A",  
  "male_gonad" = "#7E9E7E",  
  "male_tail" = "#8C7A6B"   
)

# Plot
ggplot(saturation_group, aes(x = total_reads_corrected, 
                             y = detected_transcripts,
                             group = treatment, 
                             color = treatment)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +  # use your custom fall palette
  scale_y_continuous(labels = scales::comma, limits = c(0, 75000)) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Read Depth",
    y = "Detected transcripts",
    color = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save
ggsave(filename = here::here("figures",
                             "10_satcurve_tissue.png"),
       width = 6, 
       height = 4, 
       dpi = 300)

# Plot again by proportion transcripts detected on -------------------------------
# 
# # Calculate total transcripts at fraction == 1 for each sample
# total_transcripts_per_tissue <- saturation_group[fraction == 1, .(total_transcripts = detected_transcripts), by = treatment]
# 
# # Join back and calculate proportion
# saturation_group <- merge(saturation_group, total_transcripts_per_tissue, by = "treatment")
# saturation_group[, prop_transcript := detected_transcripts / total_transcripts]
# 
# # Plot
# ggplot(saturation_group, aes(x = total_reads_corrected, 
#                              y = prop_transcript, 
#                              group = treatment, 
#                              color = treatment)) +
#   geom_line(alpha = 0.7) +
#   geom_point(size = 1.5, alpha = 0.7) +
#   scale_color_manual(values = custom_colors) +  # use your custom fall palette
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(.5, 1)) +
#   scale_x_continuous(labels = scales::comma) +
#   labs(x = "Read Depth",
#        y = "Proportion of transcripts Detected",
#        color = "Group") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# # Save
# ggsave(filename = here::here("figures",
#                              "10_satcurve_prop_tissue.png"), 
#        width = 6, 
#        height = 4, 
#        dpi = 300)

# Collapse all samples and subsample --------------------------------------

# Sum all numeric columns row-wise
num_cols <- which(sapply(count_matrix, is.numeric))
count_matrix_numeric <- as.matrix(count_matrix[, ..num_cols])
collapsed_counts <- rowSums(count_matrix_numeric, na.rm = TRUE)

# Prepare subsampling parameters
total_reads <- sum(collapsed_counts)
probs <- collapsed_counts / total_reads

# Initialize data.table to hold subsampling results
saturation_all <- data.table(
  fraction = fractions,
  total_reads = round(total_reads * fractions),
  detected_transcripts = integer(length(fractions))
)

# Perform multinomial subsampling at each fraction
saturation_all[, detected_transcripts := sapply(total_reads, function(n_reads) {
  subsampled <- rmultinom(1, size = n_reads, prob = probs)
  sum(subsampled > 0)
})]

# Correct total reads for mapping rate
saturation_all[, total_reads_corrected := round(total_reads / mapping_rate * fraction)]

# Plot 
ggplot(saturation_all, aes(x = total_reads_corrected, y = detected_transcripts)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "Read Depth",
       y = "transcripts Detected") +
  theme_minimal() +
  theme(legend.position = "right")


# Save
ggsave(filename = here::here("figures",
                             "10_satcurve_all.png"),
       width = 6, 
       height = 4, 
       dpi = 300)

# Plot again by proportion transcripts detected on --------------------------------
# 
# # Calculate proportion of transcripts detected
# saturation_all[, prop_transcript := detected_transcripts / max(detected_transcripts)]
# 
# # Plot 
# ggplot(saturation_all, aes(x = total_reads_corrected, y = prop_transcript)) +
#   geom_line(alpha = 0.7) +
#   geom_point(size = 1.5, alpha = 0.7) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(.5, 1)) +
#   scale_x_continuous(labels = scales::comma) +
#   labs(x = "Read Depth",
#        y = "Proportion of transcripts Detected") +
#   theme_minimal() +
#   theme(legend.position = "right")
# 
# # Save
# ggsave(filename = here::here("figures",
#                              "10_satcurve_prop_all.png"),
#        width = 6, 
#        height = 4, 
#        dpi = 300)
