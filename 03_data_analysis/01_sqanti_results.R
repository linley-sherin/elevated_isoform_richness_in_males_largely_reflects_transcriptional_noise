# Look at SQANTI results for isoforms

# (1) Plot # annotated gene/isoforms
# (2) look at sqanti categories
# (3) Replot after removing certain categories

library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)
library(scales)
library(colorspace)
library(patchwork)

################################################################################
### Load data
################################################################################

# Load cpm matrices for gonad and tail
cpm_gonad <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "16_edgeR_filterexp",
                                 "01_gonadtran_TMMTPMmatrix.csv"))

cpm_tail <- read.csv(here::here("data",
                                "01_build_longread_transcriptome",
                                "16_edgeR_filterexp",
                                "02_tailtran_TMMTPMmatrix.csv"))

# Function to reshape each CPM dataframe
reshape_cpm <- function(df, tissue_name) {
  df %>%
    # pivot only numeric CPM columns (e.g., X1FGP, X2MGP, etc.)
    pivot_longer(
      cols = matches("^X\\d"),
      names_to = "sample",
      values_to = "CPM"
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
    summarise(meanCPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(log2CPM = log2(meanCPM + 1))  # add log2(meanCPM + 1)
}

# Apply to both datasets
gonad_long <- reshape_cpm(cpm_gonad, "gonad")
tail_long  <- reshape_cpm(cpm_tail, "tail")

# Combine them
cpm_combined <- bind_rows(gonad_long, tail_long)

# Load SQANTI results, updated for new transcript assignments
sqanti <- read.csv(here::here("data",
                              "01_build_longread_transcriptome",
                              "16_edgeR_filterexp",
                              "SQANTI_classification_filtered.csv"))

colnames(sqanti)

# Select only columns of interest
sqanti <- sqanti %>%
  dplyr::select("gene",
         "transcript", 
         "associated_gene",
         "associated_transcript",
         "structural_category",
         "subcategory")

# Filter sqanti by genes in cpm_combined
sqanti <- sqanti %>%
  dplyr::filter(transcript %in% cpm_combined$transcript)

################################################################################
### (1) Plot # annotated gene/isoforms
################################################################################

# Gene summary
gene_summary <- sqanti %>%
  # Group by the unique ID first
  group_by(gene) %>% 
  # Pick only the first associated_gene name for each unique ID
  dplyr::slice(1) %>% 
  ungroup() %>%
  mutate(category = ifelse(grepl("novel", associated_gene, ignore.case = TRUE),
                           "novel", "annotated"),
         type = "Gene") %>%
  dplyr::count(type, category)

# Transcript summary
transcript_summary <- sqanti %>%
  # Group by the unique transcript ID
  group_by(transcript) %>% 
  # Pick the first association to avoid inflation
  dplyr::slice(1) %>% 
  ungroup() %>%
  mutate(category = ifelse(grepl("novel", associated_transcript, ignore.case = TRUE),
                           "novel", "annotated"),
         type = "Transcript") %>%
  dplyr::count(type, category)

# Combine summaries
combined_summary <- bind_rows(gene_summary, transcript_summary)

# Ensure 'annotated' is below 'novel' in stacked bars
combined_summary$category <- factor(combined_summary$category, levels = c("annotated", "novel"))

# Compute proportions and create label
combined_summary <- combined_summary %>%
  group_by(type) %>%
  mutate(prop = n / sum(n),
         label_text = paste0(round(prop * 100, 1), "% (", scales::comma(n), ")")) %>%
  ungroup()

# Plot
ggplot(combined_summary, aes(x = type, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  scale_fill_manual(values = c("annotated" = "grey70", "novel" = "#f4a1a1"))

# Save
ggsave(
  filename = here::here("figures",
                        "05_sqanti_annotated_novel.png"),
  width = 4,
  height = 4,
  dpi = 300
)

################################################################################
### (2) look at sqanti categories
################################################################################

structural_order <- c(
  "full-splice_match",
  "incomplete-splice_match",
  "novel_in_catalog",
  "novel_not_in_catalog",
  "genic",
  "genic_intron",
  "intergenic",
  "antisense",
  "fusion"
)

subcategory_order <- c(
  # green
  "reference_match",
  "alternative_3end",
  "alternative_5end",
  "alternative_3end5end",
  # blue
  "3prime_fragment",
  "5prime_fragment",
  "internal_fragment",
  # pink
  "combination_of_known_junctions",
  "combination_of_known_splicesites",
  # purple
  "at_least_one_novel_splicesite",
  # grey
  "multi-exon",
  "intron_retention"
)

subcategory_summary <- sqanti %>%
  distinct(transcript, structural_category, subcategory) %>%
  dplyr::count(structural_category, subcategory) %>%
  mutate(
    structural_category = factor(structural_category, levels = structural_order),
    subcategory = factor(subcategory, levels = subcategory_order)
  ) %>%
  group_by(structural_category) %>%
  mutate(
    prop = n / sum(n),
    label_text = paste0("n=", scales::comma(n))
  ) %>%
  ungroup()

# Base colors for structural categories
base_colors <- c(
  "incomplete-splice_match" = "#4f81bd",  # blue
  "full-splice_match"       = "#4bb14b",  # green
  "novel_in_catalog"        = "#f28e8e",  # pale red
  "novel_not_in_catalog"    = "#c27ba0"   # purple
)

# Subcategories that appear in >1 structural_category — shared subcategories
shared_subcats <- subcategory_summary %>%
  dplyr::count(subcategory) %>%
  filter(n > 1) %>%
  pull(subcategory)

# Create grey tones for shared categories
shared_colors <- setNames(
  c("grey", "lightgrey")[seq_along(shared_subcats)],
  shared_subcats
)

# Generate shades for *non-shared* subcategories
subcategory_colors <- subcategory_summary %>%
  filter(!subcategory %in% shared_subcats) %>%
  group_by(structural_category) %>%
  summarise(subcats = list(unique(subcategory)), .groups = "drop") %>%
  rowwise() %>%
  mutate(colors = list(
    lighten(base_colors[structural_category],
            seq(0, 0.4, length.out = length(subcats)))
  )) %>%
  unnest(c(subcats, colors)) %>%
  dplyr::select(subcategory = subcats, color = colors)

# Combine structural-category–based colors with greys
final_colors <- c(
  setNames(subcategory_colors$color, subcategory_colors$subcategory),
  shared_colors
)

subcategory_colors$subcategory <- factor(subcategory_colors$subcategory,
                                         levels = subcategory_order)
# Plot
ggplot(subcategory_summary, aes(x = structural_category, y = n, fill = subcategory)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), color = "black", size = 3.5) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = final_colors) +
  labs(x = NULL,
       y = "Transcript count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save
ggsave(
  filename = here::here("figures",
                        "05_sqanti_structural_categories.png"),
  width = 9,
  height = 6,
  dpi = 300
)

################################################################################
### (3) Replot after removing certain categories
################################################################################

# Set plot theme (font size = 8)
theme_set(theme_bw(base_size = 8))

glimpse(sqanti)

# Remove all isoforms classified as genic, genic_intron, intergenic, antisense, fusion
sqanti_filtered <- sqanti %>%
  filter(!structural_category %in% c("genic", "genic_intron", "intergenic", "antisense", "fusion"))

# Gene summary
gene_summary <- sqanti_filtered %>%
  # Group by the unique ID first
  group_by(gene) %>% 
  # Pick only the first associated_gene name for each unique ID
  dplyr::slice(1) %>% 
  ungroup() %>%
  mutate(category = ifelse(grepl("novel", associated_gene, ignore.case = TRUE),
                           "novel", "annotated"),
         type = "Gene") %>%
  dplyr::count(type, category)

# Transcript summary
transcript_summary <- sqanti_filtered %>%
  # Group by the unique transcript ID
  group_by(transcript) %>% 
  # Pick the first association to avoid inflation
  dplyr::slice(1) %>% 
  ungroup() %>%
  mutate(category = ifelse(grepl("novel", associated_transcript, ignore.case = TRUE),
                           "novel", "annotated"),
         type = "Transcript") %>%
  dplyr::count(type, category)

# Combine summaries
combined_summary <- bind_rows(gene_summary, transcript_summary)

# Ensure 'annotated' is below 'novel' in stacked bars
combined_summary$category <- factor(combined_summary$category, levels = c("annotated", "novel"))

# Compute proportions and create label
combined_summary <- combined_summary %>%
  group_by(type) %>%
  mutate(prop = n / sum(n),
         label_text = paste0(round(prop * 100, 1), "% (", scales::comma(n), ")")) %>%
  ungroup()

# Plot
novelty_plot <- ggplot(combined_summary, aes(x = type, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  scale_fill_manual(values = c("annotated" = "grey70", "novel" = "#f4a1a1"))

structural_order <- c(
  "full-splice_match",
  "incomplete-splice_match",
  "novel_in_catalog",
  "novel_not_in_catalog"
)

subcategory_order <- c(
  # green
  "reference_match",
  "alternative_3end",
  "alternative_5end",
  "alternative_3end5end",
  # blue
  "3prime_fragment",
  "5prime_fragment",
  "internal_fragment",
  # pink
  "combination_of_known_junctions",
  "combination_of_known_splicesites",
  # purple
  "at_least_one_novel_splicesite",
  # grey
  "multi-exon",
  "intron_retention"
)

subcategory_summary <- sqanti_filtered %>%
  distinct(transcript, structural_category, subcategory) %>%
  dplyr::count(structural_category, subcategory) %>%
  mutate(
    structural_category = factor(structural_category, levels = structural_order),
    subcategory = factor(subcategory, levels = subcategory_order)
  ) %>%
  group_by(structural_category) %>%
  mutate(
    prop = n / sum(n),
    label_text = paste0("n=", scales::comma(n))
  ) %>%
  ungroup()

# Base colors for structural categories
base_colors <- c(
  "incomplete-splice_match" = "#4f81bd",  # blue
  "full-splice_match"       = "#4bb14b",  # green
  "novel_in_catalog"        = "#f28e8e",  # pale red
  "novel_not_in_catalog"    = "#c27ba0"   # purple
)

# Subcategories that appear in >1 structural_category — shared subcategories
shared_subcats <- subcategory_summary %>%
  dplyr::count(subcategory) %>%
  filter(n > 1) %>%
  pull(subcategory)

# Create grey tones for shared categories
shared_colors <- setNames(
  c("grey", "lightgrey")[seq_along(shared_subcats)],
  shared_subcats
)

# Generate shades for *non-shared* subcategories
subcategory_colors <- subcategory_summary %>%
  filter(!subcategory %in% shared_subcats) %>%
  group_by(structural_category) %>%
  summarise(subcats = list(unique(subcategory)), .groups = "drop") %>%
  rowwise() %>%
  mutate(colors = list(
    lighten(base_colors[structural_category],
            seq(0, 0.4, length.out = length(subcats)))
  )) %>%
  unnest(c(subcats, colors)) %>%
  dplyr::select(subcategory = subcats, color = colors)

# Combine structural-category–based colors with greys
final_colors <- c(
  setNames(subcategory_colors$color, subcategory_colors$subcategory),
  shared_colors
)

subcategory_colors$subcategory <- factor(subcategory_colors$subcategory,
                                         levels = subcategory_order)
# Plot
cat_plot <- ggplot(subcategory_summary, aes(x = structural_category, y = n, fill = subcategory)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5), color = "black", size = 3) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = final_colors) +
  labs(x = NULL,
       y = "Transcript count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

patched <- (novelty_plot | cat_plot) + 
  plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))

patched

# Save
ggsave(plot = patched,
       here::here("figures",
                  "05_combined_sqanti_plots.png"),
  width = 14,
  height = 9,
  dpi = 300
)

################################################################################
### (4) UPDATE GTF, TPM and COUNT MATRICES SO MATCH sqanti_w_orf
################################################################################
# 
# # Load TPM matrices for gonad and tail
# TPM_gonad <- read.csv(here::here("data",
#                                  "isoseq",
#                                  "10_edgeR_filterexp",
#                                  "01_gonadtran_TMMTPMmatrix.csv")) %>%
#   filter(transcript %in% sqanti_filtered$transcript)
# TPM_tail <- read.csv(here::here("data",
#                                 "isoseq",
#                                 "10_edgeR_filterexp",
#                                 "02_tailtran_TMMTPMmatrix.csv")) %>%
#   filter(transcript %in% sqanti_filtered$transcript)
# count_gonad <- read.csv(here::here("data",
#                                    "isoseq",
#                                    "10_edgeR_filterexp",
#                                    "01_gonadtran_countmatrix.csv")) %>%
#   filter(transcript %in% sqanti_filtered$transcript)
# count_tail <- read.csv(here::here("data",
#                                   "isoseq",
#                                   "10_edgeR_filterexp",
#                                   "02_tailtran_countmatrix.csv")) %>%
#   filter(transcript %in% sqanti_filtered$transcript)
# 
# # Resave
# write.csv(TPM_gonad,
#           here::here("data","01_build_longread_transcriptome","16_edgeR_filterexp",
#                      "01_gonadtran_TMMTPMmatrix_long_sqantifiltered.csv"),
#           row.names = FALSE)
# write.csv(TPM_tail,
#           here::here("data","01_build_longread_transcriptome","16_edgeR_filterexp",
#                      "02_tailtran_TMMTPMmatrix_long_sqantifiltered.csv"),
#           row.names = FALSE)
# write.csv(count_gonad,
#           here::here("data","01_build_longread_transcriptome","16_edgeR_filterexp",
#                      "01_gonadtran_countmatrix_long_sqantifiltered.csv"),
#           row.names = FALSE)
# write.csv(count_tail,
#           here::here("data","01_build_longread_transcriptome","16_edgeR_filterexp",
#                      "02_tailtran_countmatrix_long_sqantifiltered.csv"),
#           row.names = FALSE)
# 
# # Load gtfs
# gtf <- import(here::here("data",
#                          "GTF",
#                          "04_filtered_by_exp_long.gtf"))
# 
# # Filter gtf
# gtf_sqanti <- gtf[gtf$transcript_id %in% sqanti_filtered$transcript]
# 
# # Save
# rtracklayer::export(
#   gtf_sqanti,
#   here::here("data",
#              "GTF",
#              "05_filtered_by_exp_long_sqantifiltered.gtf")
# )
