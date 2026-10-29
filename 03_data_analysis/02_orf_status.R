# Investigate ORF status (whether isoform has a reading frame) and how
# isoform features relate to expression

### FILTER NOTICE: I REMOVE ALL GENES THAT HAVE ONLY 1 ISOFORM IN A TISSUE
### (ONLY FOR THE FIGURE)

# (1) Determine ORF status of all transcripts
# (2) Plot as in benetiere (against rel abundance)
# (3) Plot ORF_status vs. isoform absolute abundance
# (4) Plot number of isoforms per gene, male vs. female (scatterplot)

library(GenomicFeatures)
library(Biostrings)
library(dplyr)
library(stringr)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(PNWColors)
library(purrr)
library(tidyverse)
library(here)
library(gghalves)
library(ggbeeswarm)
library(forcats)
library(broom)

################################################################################
### (1) Determine ORF status of all transcripts
################################################################################

### Load data

# Load TPM matrices for gonad and tail
TPM_gonad <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "16_edgeR_filterexp",
                                 "01_gonadtran_TMMTPMmatrix.csv"))

TPM_tail <- read.csv(here::here("data",
                                "01_build_longread_transcriptome",
                                "16_edgeR_filterexp",
                                "02_tailtran_TMMTPMmatrix.csv"))

# Function to reshape each TPM dataframe
reshape_TPM <- function(df, tissue_name) {
  df %>%
    # pivot only numeric TPM columns (e.g., X1FGP, X2MGP, etc.)
    pivot_longer(
      cols = matches("^X\\d"),
      names_to = "sample",
      values_to = "TPM"
    ) %>%
    # Extract sex info from column names (F = female, M = male)
    mutate(
      sex = case_when(
        str_detect(sample, "F") ~ "female",
        str_detect(sample, "M") ~ "male",
        TRUE ~ NA_character_
      ),
      tissue = tissue_name
    ) %>%
    group_by(transcript, gene, tissue, sex) %>% # Group by each treatment identify and calculate mean
    summarise(mean_isoformTPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(log2_isoformTPM = log2(mean_isoformTPM + 1))  # add log2(mean_isoformTPM + 1)
}

# Apply to both datasets
gonad_long <- reshape_TPM(TPM_gonad, "gonad")
tail_long  <- reshape_TPM(TPM_tail, "tail")

# Combine them
TPM_combined <- bind_rows(gonad_long, tail_long)

# Collapse TPM within each treatment group (take avg TPM including 0s)
# Find relative abundance within each treatment
TPM_combined <- TPM_combined %>%
  group_by(gene, tissue, sex) %>%
  mutate(
    geneTPM = sum(mean_isoformTPM, na.rm = TRUE),
    relAbundance = mean_isoformTPM / geneTPM
  ) %>%
  ungroup()

# Load output from sqanti
class <- read.csv(here::here("data",
                             "01_build_longread_transcriptome",
                             "16_edgeR_filterexp",
                             "SQANTI_classification_filtered.csv"),
                  header = T)

class <- class %>%
  filter(transcript %in% TPM_combined$transcript)

# Load ORF output from transdecoder
ORF <- read.csv(here::here("data",
                           "01_build_longread_transcriptome",
                           "13_transdecoder",
                           "parsed_results.csv"))

# Load gene reassignment map
map <- read.csv(here::here("data",
                           "01_build_longread_transcriptome",
                           "14_remove_redundancy",
                           "rosetta_stone.csv"))

# Join map with ORF
ORF <- ORF %>%
  left_join(map, by = c("transcript_id", "gene_id")) %>%
  dplyr::select(-transcript_id, -gene_id, -orf_id) %>%
  dplyr::rename("transcript_id" = "new_transcript_id", 
                "gene_id" = "new_gene_id")

# Join with class
class_with_orf <- class %>%
  left_join(ORF, by = c("transcript" = "transcript_id", "gene" = "gene_id"))

### Determine if isoform has evidence of NMD

# Evidence of NMD = a stop codon >50bp upstream of the last exon-exon boundary

# Load genome-aligned GTF with exon structures (you could also use gff3 from transdecoder workflow)
gtf <- import(here::here("data",
                         "GTF",
                         "04_filtered_by_exp_long.gtf"))
gtf_df <- as.data.frame(gtf)

# Filter by TPM_all transcripts
gtf_df <- gtf_df %>%
  filter(transcript_id %in% TPM_combined$transcript)

# Filter for exons
exons <- gtf_df %>% 
  filter(type == "exon") %>%
  dplyr::select(seqnames, start, end, strand, gene_id, transcript_id)

# Get list of transcripts that appear in both exons and ORF_filtered
tx_ids <- intersect(unique(exons$transcript_id), ORF$transcript_id)

##### Find the last junction for each transcript (depends on strand)

# Split into list of data.frames, one per transcript
split_exons <- exons %>%
  group_by(transcript_id) %>%
  group_split()

# Function to get penult junction for one transcript
get_last_junction <- function(df) {
  strand <- as.character(df$strand[1])
  n_exons <- nrow(df)
  if (n_exons < 2) return(NULL)
  
  if (strand == "+") {
    df_ordered <- df %>% arrange(start)
    penult_end <- df_ordered$end[n_exons - 1]
  } else if (strand == "-") {
    df_ordered <- df %>% arrange(desc(start))
    penult_end <- df_ordered$start[2]
  } else {
    return(NULL)
  }
  
  return(data.frame(transcript_id = df$transcript_id[1],
                    strand = strand,
                    penult_exon_end = penult_end))
}

# Apply to all transcripts
last_junctions <- map_dfr(split_exons, get_last_junction)

# Join with ORF CDS ends
cds_nmd <- ORF %>%
  dplyr::select(transcript_id, cds_start, cds_end) %>%
  right_join(last_junctions, by = c("transcript_id")) %>%
  mutate(
    stop_codon_pos = case_when(
      strand == "+" ~ cds_end,
      strand == "-" ~ cds_start,
      TRUE ~ NA_real_
    ),
    dist_to_junction = case_when(
      strand == "+" ~ penult_exon_end - stop_codon_pos,
      strand == "-" ~ stop_codon_pos - penult_exon_end,
      TRUE ~ NA_real_
    ),
    nmd_candidate = dist_to_junction > 50
  )

# Check results
table(cds_nmd$nmd_candidate, cds_nmd$strand)

# Save output
write_csv(cds_nmd, here::here("data",
                              "01_build_longread_transcriptome",
                              "13_transdecoder",
                              "nmd_prediction_results.csv"))

#### REMINDER: this will 'predict' NMD for all orf_types, but only complete and 5prime_partials can actually be NMD
#### internal or 3prime_partials don't have stop codons and so it's not actually possible to determine NMD

# Join with class_with_orf
class_with_orf_nmd <- class_with_orf %>%
  dplyr::select(-cds_start, -cds_end) %>%
  left_join(cds_nmd, by = c("transcript" = "transcript_id"))

### Designate ORF_status

# Assign orf metrics
class_with_orf_nmd <- class_with_orf_nmd %>%
  mutate(
    orf_status = case_when(
      subcategory == "intron_retention" ~ "intron_retention", # First true condition wins
      orf_type == "complete" & nmd_candidate == F ~ "full_orf",
      orf_type == "complete" & nmd_candidate == T ~ "NMD_candidate",
      orf_type == "5prime_partial" & nmd_candidate == F ~ "no_orf",
      orf_type == "5prime_partial" & nmd_candidate == T ~ "no_orf",
      orf_type == "3prime_partial" ~ "no_orf",
      orf_type == "internal" ~ "no_orf",
      is.na(orf_type) ~ "no_orf"),
    ptc_status = case_when(
      orf_status == "full_orf" ~ "coding",
      orf_status == "NMD_candidate" ~ "nmd",
      orf_status == "no_orf" ~ "non-coding",
      orf_status == "intron_retention" ~ "intron_retention"),
    cds_identity = factor(
      cds_seq,  # Do not subset cds_seq, keep the whole column
      levels = unique(cds_seq[!is.na(cds_seq)]),  # Only consider non-NA sequences for levels
      labels = paste("cds", seq_along(unique(cds_seq[!is.na(cds_seq)])), sep = "_")),  # Labels for unique non-NA sequences
    cds_identity = if_else(ptc_status == "coding", cds_identity, NA_character_), # Remove cds_identity for non-coding transcripts
    orf_identity = factor(
      orf_seq,  # Do not subset orf_seq, keep the whole column
      levels = unique(orf_seq[!is.na(orf_seq)]),  # Only consider non-NA sequences for levels
      labels = paste("orf", seq_along(unique(orf_seq[!is.na(orf_seq)])), sep = "_")),  # Labels for unique non-NA sequences
    orf_identity = if_else(ptc_status == "coding", orf_identity, NA_character_) # Remove orf_identity for non-coding transcripts
  )

# Sanity check: Check if any orf_id appears in more than one gene
# sanity_check <- class_with_orf_nmd %>%
#   group_by(orf_identity) %>%
#   summarise(gene_count = n_distinct(gene))

orf_status_tosave <- class_with_orf_nmd %>%
  dplyr::select(gene, transcript, orf_status, ptc_status, cds_identity, orf_identity)

write.csv(orf_status_tosave,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "orf_status.csv"))

# Clean up class_with_orf_nmd
colnames(class_with_orf_nmd)

class_with_orf_nmd <- class_with_orf_nmd %>%
  dplyr::select(gene,
                transcript, 
                chrom,
                length,
                exons, 
                structural_category,
                associated_gene,
                associated_transcript,
                ref_length,
                ref_exons,
                FSM_class,
                orf_type,
                orf_len,
                orf_seq,
                cds_start,
                cds_end,
                nmd_candidate,
                orf_status,
                ptc_status,
                cds_identity,
                orf_identity)

write.csv(class_with_orf_nmd,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "16_edgeR_filterexp",
                     "sqanti_w_orf.csv"))

################################################################################
### (2) Plot as in Benetiere et al. (against rel abundance) - OLD WAY
### Using collapsed TPM matrix where all replicates are combined within 
### a condition. 
### No additional filtering on the really lowly expressed isoforms
### across all treatments
################################################################################

# Join with orf_status and remove all transcripts with TPM = 0
TPM_combined <- TPM_combined %>%
  right_join(
    dplyr::select(class_with_orf_nmd, transcript, gene, orf_status),
    by = c("transcript", "gene")
  ) %>%
  filter(mean_isoformTPM > 0) %>%
  mutate(
    tissue = factor(tissue, levels = c("tail", "gonad"))
  )

### FILTER NOTICE: I REMOVE ALL GENES THAT HAVE ONLY 1 ISOFORM IN A TISSUE
# Remove transcripts with total relAbundance = 1 across all samples within a tissue (single isoform genes)
TPM_combined <- TPM_combined %>%
    group_by(transcript, tissue) %>%
    filter(all(relAbundance == 1) == FALSE) %>%
    ungroup()

# Bin relAbundance into 100 fixed-width bins (1% each)
TPM_binned <- TPM_combined %>%
  mutate(
    bin = cut(
      relAbundance,
      breaks = seq(0, 1, by = 0.05),
      include.lowest = TRUE,
      right = FALSE
    )
  )

# Summarize ORF-status proportions within each bin, tissue, and sex
bin_summary <- TPM_binned %>%
  group_by(tissue, sex, bin, orf_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sex, bin) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Clean up bin labels (optional, nicer on x-axis)
bin_summary <- bin_summary %>%
  mutate(
    bin_label = as.character(bin),
    bin_mid = as.numeric(sub("\\[(.*),.*", "\\1", bin_label)) + 0.005  # midpoint for plotting
  )

# Order ORF categories for consistent coloring
bin_summary$orf_status <- factor(
  bin_summary$orf_status,
  levels = rev(c(
    "full_orf",
    "NMD_candidate",
    "intron_retention",
    "no_orf"
  ))
)

bin_summary <- bin_summary %>%
  mutate(
    tissue = factor(tissue, levels = c("tail", "gonad"))
  )

# Plot — fixed bins (5% increments), faceted by tissue × sex
### Note, this gives error because it's not plotting isoforms that has 0 expression
ggplot(bin_summary, aes(x = bin_mid, y = proportion, fill = orf_status)) +
  geom_bar(stat = "identity", position = "fill", width = 0.05) +
  scale_x_continuous(
    name = "Isoform relative abundance",
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  scale_fill_manual(values = pnw_palette("Sailboat", 4), name = "ORF Status") +
  labs(y = "Proportion of isoforms") +
  facet_grid(sex ~ tissue) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

# Save
ggsave(
  filename = here::here("figures",
                        "06_orfstatus_relative_facetedbytreatment.png"),
  width = 5.5,           # width in inches
  height = 4,          # height in inches
  dpi = 300,           # resolution (use 300 for publication-quality)
)

################################################################################
### (3) Plot ORF_status vs. isoform absolute abundance 
################################################################################

# Define number bins
num_bins <- 20

# Bin log2_isoformTPM using equal-frequency bins within tissue and sex
TPM_binned_log <- TPM_binned %>%
  filter(mean_isoformTPM > 0) %>% # redundant but just another check
  group_by(tissue, sex) %>%
  mutate(bin = ntile(log2_isoformTPM, num_bins)) %>%
  ungroup()

# Compute median log2_isoformTPM per bin for labels
bin_medians <- TPM_binned_log %>%
  group_by(tissue, sex, bin) %>%
  summarise(bin_mid = median(log2_isoformTPM, na.rm = TRUE), .groups = "drop")

# Summarize ORF-status proportions per bin
bin_summary_log <- TPM_binned_log %>%
  group_by(tissue, sex, bin, orf_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue, sex, bin) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  left_join(bin_medians, by = c("tissue", "sex", "bin"))

# Order ORF categories
bin_summary_log$orf_status <- factor(
  bin_summary_log$orf_status,
  levels = rev(c(
    "full_orf",
    "NMD_candidate",
    "intron_retention",
    "no_orf"
  )))

bin_summary_log <- bin_summary_log %>%
  mutate(
    tissue = factor(tissue, levels = c("tail", "gonad"))
  )

# Plot with uniform bin widths, x-axis labeled by median log2_isoformTPM
ggplot(bin_summary_log, aes(x = factor(bin), y = proportion, fill = orf_status)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  scale_x_discrete(
    name = "Isoform expression (log2 TPM)",
    labels = function(b) {
      # Convert factor levels to numeric
      b_num <- as.numeric(b)
      # For each tissue/sex combination, take median bin value
      # (for simplicity here, just use overall median across all tissue/sex)
      sapply(b_num, function(x) round(median(bin_medians$bin_mid[bin_medians$bin == x]), 2))
    }
  ) +
  scale_fill_manual(values = pnw_palette("Sailboat", 4), name = "ORF Status") +
  labs(y = "Proportion of isoforms") +
  facet_grid(sex ~ tissue) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(face = "bold"))

# Save figure
ggsave(
  filename = here::here("figures", 
                        "06_orfstatus_absolute_facetedbytreatment.png"),
  width = 5.5,
  height = 4,
  dpi = 300
)

################################################################################
### (4) Plot number of isoforms per gene, male vs. female (scatterplot)
################################################################################

# Count number of expressed isoforms per gene, separately by tissue and sex
isoform_counts <- TPM_combined %>%
  filter(mean_isoformTPM > 0) %>%  # only expressed isoforms
  group_by(tissue, sex, gene) %>%
  summarise(num_isoforms = n(), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = num_isoforms, values_fill = 0)

# Set tissue order: tail first, gonad second
isoform_counts <- isoform_counts %>%
  mutate(tissue = factor(tissue, levels = c("tail", "gonad")))

# Find the maximum isoform count to set axes
max_isoforms <- max(isoform_counts$female, isoform_counts$male, na.rm = TRUE)

# Scatterplot with fixed axes
ggplot(isoform_counts, aes(x = female, y = male)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ tissue, scales = "fixed") +  # fixed scales
  scale_x_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
  labs(
    x = "Female Isoform Count",
    y = "Male Isoform Count",
  ) +
  theme_minimal(base_size = 12)

ggplot(isoform_counts, aes(x = female, y = male)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # trend line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ tissue, scales = "fixed") +
  scale_x_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
  labs(
    x = "Female Isoform Count",
    y = "Male Isoform Count"
  ) +
  theme_minimal(base_size = 12)

# Save figure
ggsave(
  filename = here::here("figures", 
                        "06_num_isoforms_malevfemale.png"),
  width = 6,
  height = 4,
  dpi = 300
)

# ### Remove genes with only 1 isoform within a tissue
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(scales)
# 
# # 1. Count and Filter
# isoform_counts2 <- TPM_combined %>%
#   filter(mean_isoformTPM > 0) %>%
#   group_by(tissue, sex, gene) %>%
#   summarise(num_isoforms = n(), .groups = "drop") %>%
#   # Ensure genes expressed in only one sex are kept as 0 in the other
#   complete(gene, nesting(tissue, sex), fill = list(num_isoforms = 0)) %>%
#   # Filter: Keep only genes that have > 1 isoform TOTAL in that tissue (sum of both sexes)
#   group_by(tissue, gene) %>%
#   filter(sum(num_isoforms) > 1) %>% 
#   ungroup() %>%
#   # Pivot for plotting
#   pivot_wider(names_from = sex, values_from = num_isoforms, values_fill = 0)
# 
# # 2. Set tissue order
# isoform_counts2 <- isoform_counts2 %>%
#   mutate(tissue = factor(tissue, levels = c("tail", "gonad")))
# 
# # 4. Plot
# ggplot(isoform_counts2, aes(x = female, y = male)) +
#   geom_point(alpha = 0.2, size = 0.5) +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") +  # trend line
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   facet_wrap(~ tissue, scales = "fixed") +
#   scale_x_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, max_isoforms + 1), expand = c(0, 0)) +
#   labs(
#     x = "Female Isoform Count",
#     y = "Male Isoform Count"
#   ) +
#   theme_minimal(base_size = 12)
