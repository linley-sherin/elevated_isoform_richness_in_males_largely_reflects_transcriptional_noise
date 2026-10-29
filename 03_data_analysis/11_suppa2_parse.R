# Parse results of suppa2

# (1) Load TPM
# (2) Match transcripts and venn diagram
# (3) Check whether transcripts Short Only-reads have unique junctions
# (4) Check whether transcript TPM in each dataset
# (5) Load and compare events discovered
# (6) Parse SUPPA2 results
# (7) Plot all dpsi before coloured by sig
# (8) Stacked bar plots of events 

library(multcompView)
library(emmeans)
library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)
library(ggrepel)
library(GenomicRanges)
library(ggvenn)
library(grid)
library(fs)
library(forcats)
library(stringr)
library(purrr)
library(patchwork)

################################################################################
### (1) Load TPM
################################################################################

source("load_mt.R")

# Load tpm
tpm_short_tail <- read.csv(here::here("data", 
                                      "02_build_shortread_transcriptome", 
                                      "05_edgeR",
                                      "02_tailtran_TMMTPMmatrix_short.csv"))[,-1]
tpm_short_gonad <- read.csv(here::here("data", 
                                       "02_build_shortread_transcriptome", 
                                       "05_edgeR",
                                       "01_gonadtran_TMMTPMmatrix_short.csv"))[,-1]
tpm_long_tail <- read.csv(here::here("data", 
                                     "01_build_longread_transcriptome", 
                                     "16_edgeR_filterexp",
                                     "02_tailtran_TMMTPMmatrix.csv"))[,-1]
tpm_long_gonad <- read.csv(here::here("data", 
                                      "01_build_longread_transcriptome", 
                                      "16_edgeR_filterexp",
                                      "01_gonadtran_TMMTPMmatrix.csv"))


# define a function to transform each dataset
process_tpm <- function(df, type_val, tissue_val) {
  df %>%
    # Convert from wide to long: move all sample columns (X1FTP, etc.) into two columns
    pivot_longer(
      cols = starts_with("X"), 
      names_to = "sample_id", 
      values_to = "tpm"
    ) %>%
    mutate(
      data_type = type_val,
      tissue = tissue_val,
      # Parse the sample_id (e.g., "X1FTP")
      sex = ifelse(substr(sample_id, 3, 3) == "F", "female", "male")
    ) %>%
    # Clean up: remove the raw sample_id if no longer needed
    dplyr::select(transcript, gene, data_type, tissue, sex, sample_id, tpm)
}

# apply the function to each of your 4 dataframes
short_tail  <- process_tpm(tpm_short_tail, "short", "tail")
long_gonad  <- process_tpm(tpm_long_gonad, "long", "gonad")
short_gonad <- process_tpm(tpm_short_gonad, "short", "gonad")
long_tail   <- process_tpm(tpm_long_tail, "long", "tail")

# combine everything into one master long-format dataframe
tpm_all_long <- bind_rows(
  short_tail, 
  long_gonad, 
  short_gonad, 
  long_tail
)

# Clean up
rm(short_gonad)
rm(short_tail)
rm(long_gonad)
rm(long_tail)
rm(process_tpm)
rm(mt_all)
rm(mt_tissue)
rm(mt_treat)

################################################################################
### (2) Match transcripts and venn diagram
################################################################################

# Genes are already matched by reference id
# Transcripts have been matched with fuzzy_venn.py but the output
# needs to be processed

# long-read gene to transcript lookup table
lr_lookup <- bind_rows(tpm_long_tail, tpm_long_gonad) %>%
  distinct(transcript, gene) %>%
  dplyr::rename(LR_tran = transcript, LR_gene = gene)

length(unique(lr_lookup$LR_gene))
# 13,721
length(unique(lr_lookup$LR_tran))
# 64,074

# short-read gene to transcript lookup table
sr_lookup <- bind_rows(tpm_short_tail, tpm_short_gonad) %>%
  distinct(transcript, gene) %>%
  dplyr::rename(SR_tran = transcript, SR_gene = gene)

length(unique(sr_lookup$SR_gene))
# 13,945
length(unique(sr_lookup$SR_tran))
# 24,532

# Genes found in BOTH datasets
length(unique(lr_lookup$LR_gene))
length(unique(sr_lookup$SR_gene))
shared_genes <- intersect(lr_lookup$LR_gene, sr_lookup$SR_gene) %>%
  unique()

# Genes found ONLY in Long-Read (missing from Short-Read)
unique_to_lr <- setdiff(lr_lookup$LR_gene, sr_lookup$SR_gene)

# Genes found ONLY in Short-Read (missing from Long-Read)
unique_to_sr <- setdiff(sr_lookup$SR_gene, lr_lookup$LR_gene)

# Quick Summary
cat("Summary Table:",
    "\nShared:       ", length(shared_genes),
    "\nUnique to LR: ", length(unique_to_lr),
    "\nUnique to SR: ", length(unique_to_sr))

# Read in matched csv and filter for only shared genes
matched_tran <- read.csv(here::here("data",
                                    "03_data_analysis",
                                    "10_match_transcripts",
                                    "0_matched_ids.csv")) %>%
  filter(gene_id %in% shared_genes)

# Check status
table(matched_tran$status)

### Replicate bar plot for events (below)

# Summarize the raw isoform statuses from your matched_tran table
plot_data_iso_match <- matched_tran %>%
  group_by(status) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(
    match_status = case_when(
      status == "lr_unique" ~ "Long",
      status == "matched"   ~ "Shared",
      status == "sr_unique" ~ "Short"
    ),
    match_status = factor(match_status, levels = c("Long", "Shared", "Short"))
  )

total_isoforms <- sum(plot_data_iso_match$count, na.rm = TRUE)
y_max_iso <- max(plot_data_iso_match$count, na.rm = TRUE) * 1.25

p_iso_concordance <- ggplot(plot_data_iso_match, aes(x = match_status, y = count)) +
  # Set fill directly to a single gray color instead of mapping it to a variable
  geom_bar(stat = "identity", width = 0.85, fill = "grey60", color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = paste0(round((count / total_isoforms) * 100, 1), "%")),
    vjust = -0.5,       
    size = 3.5,         
    color = "#2c3e50"   
  ) +
  scale_y_continuous(labels = scales::comma, limits = c(0, y_max_iso), expand = c(0, 0)) +
  labs(
    y = NULL,  
    x = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(), 
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
    legend.position = "none" 
  )

p_iso_concordance

# clean up
rm(lr_lookup)
rm(sr_lookup)
rm(tpm_long_gonad)
rm(tpm_long_tail)
rm(tpm_short_gonad)
rm(tpm_short_tail)
rm(unique_to_lr)
rm(unique_to_sr)

################################################################################
### (3) Check whether transcripts with short-only-reads have unique junctions
################################################################################

# Use df of matched trans, already filtered for shared_genes
glimpse(matched_tran)

# Load the GTF files
LRT <- rtracklayer::import(here::here("data",
                                      "GTF",
                                      "04_filtered_by_exp_long.gtf")) %>%
  data.frame() %>%
  filter(gene_id %in% matched_tran$gene_id)

SRT <- rtracklayer::import(here::here("data",
                                      "GTF",
                                      "06_shortreadpipe_filtered_by_exp_short.gtf")) %>%
  data.frame() %>%
  filter(gene_id %in% matched_tran$gene_id)

# Sanity check
cat("LRT contains", length(unique(LRT$gene_id)), "genes.\n")
cat("LRT contains", length(unique(LRT$transcript_id)), "isoforms.\n")
cat("LRT contains", length(unique(LRT$gene_id)), "genes.\n")
cat("SRT contains", length(unique(SRT$transcript_id)), "isoforms.\n")

glimpse(matched_tran)
glimpse(LRT)
glimpse(SRT)

# Get junctions per transcript
get_tx_junctions <- function(df) {
  df %>%
    filter(type == "exon") %>%
    arrange(transcript_id, start) %>%
    group_by(transcript_id, gene_id) %>%
    mutate(
      next_start = lead(start),
      donor = end,
      acceptor = next_start
    ) %>%
    filter(!is.na(next_start)) %>%
    transmute(
      transcript_id,
      gene_id,
      junction = paste(seqnames, donor, acceptor, strand, sep = ":")
    )
}

lrt_junc <- get_tx_junctions(LRT)
srt_junc <- get_tx_junctions(SRT)

# Create global junction set for LRT 
lrt_junc_set <- unique(lrt_junc$junction)

# Look only at sr_unique isoforms
srt_unique_tx <- matched_tran %>%
  filter(status == "sr_unique") %>%
  pull(short_read_id) %>%
  unique()

# Classify each SRT sr_unique isoform
srt_classification <- srt_junc %>%
  filter(transcript_id %in% srt_unique_tx) %>%
  mutate(junc_in_lrt = junction %in% lrt_junc_set) %>%
  group_by(transcript_id) %>%
  summarise(
    total_junc = n(),
    novel_junc = sum(!junc_in_lrt),
    has_novel_junction = novel_junc > 0,
    .groups = "drop"
  )

# Check results
srt_classification %>%
  dplyr::count(has_novel_junction)

################################################################################
### (4) Check transcript TPM in each dataset
################################################################################

set_theme(theme_classic() + theme(strip.background = element_rect(fill = 'grey90', color = NA)))

# Goal:
# Create a scatterplot comparing long-read (x-axis) vs short-read (y-axis)
# transcript expression, and colour points by transcript class:
# - matched
# - lr-unique
# - sr-unique (novel junctions)
# - sr-unique (novel assembly)

# Collapse TPM per transcript
tpm_wide <- tpm_all_long %>%
  group_by(transcript, data_type) %>%
  summarise(mean_tpm = mean(tpm, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = data_type,
    values_from = mean_tpm
  )

# Build transcript-level bridge using matched_tran
tx_map <- matched_tran %>%
  dplyr::select(short_read_id, long_read_id, status)

# Attach SRT TPM
srt_tpm <- tpm_wide %>%
  transmute(
    short_read_id = transcript,
    short_tpm = short
  )

# Attach LRT TPM
lrt_tpm <- tpm_wide %>%
  transmute(
    long_read_id = transcript,
    long_tpm = long
  )

### Prepare to make a scatterplot

# Combine everything
plot_df <- tx_map %>%
  left_join(srt_tpm, by = "short_read_id") %>%
  left_join(lrt_tpm, by = "long_read_id")

# Add SR isoform structure class
plot_df <- plot_df %>%
  left_join(
    srt_classification %>%
      dplyr::select(transcript_id, has_novel_junction),
    by = c("short_read_id" = "transcript_id")
  ) %>%
  mutate(
    class = case_when(
      status == "matched" ~ "Shared",
      status == "lr_unique" ~ "Long Only",
      status == "sr_unique" & has_novel_junction ~ "Novel Junction",
      status == "sr_unique" & !has_novel_junction ~ "Novel Assembly",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(class))

# Handle missing values for plotting
plot_df <- plot_df %>%
  mutate(
    short_tpm = ifelse(is.na(short_tpm), 0, short_tpm),
    long_tpm  = ifelse(is.na(long_tpm), 0, long_tpm)
  )

my_palette <- c(
  "Shared" = "#999999",
  "Novel Assembly" = "#1B9E77",
  "Novel Junction" = "#D95F02"#,
  #"Long Only" = "#7570B3"
)

breaks_shared <- c(0, 1, 5, 20, 100, 1000, 5000)

# Collapse SR TPM per transcript
sr_expr <- tpm_all_long %>%
  filter(data_type == "short") %>%
  group_by(transcript) %>%
  summarise(sr_tpm = mean(tpm, na.rm = TRUE), .groups = "drop")

# Build plotting table
box_df <- plot_df %>%
  filter(class != "Long Only") %>%
  mutate(class = factor(class, levels = c(
    "Shared", 
    "Novel Assembly", # Green - Now in the middle
    "Novel Junction" # Orange - Now on the right
  )))

# log transform for model
box_df <- box_df %>%
  mutate(log_tpm = log1p(short_tpm))

# # ANOVA model
# model <- aov(log_tpm ~ class, data = box_df)
# 
# # Estimated marginal means
# emm <- emmeans(model, ~ class)
# 
# # Compact letter display
# cld_results <- cld(
#   emm,
#   Letters = letters,
#   adjust = "tukey"
# )
# 
# cld_df <- cld_results %>%
#   as.data.frame() %>%
#   mutate(
#     .group = case_when(
#       class == "Shared" ~ "a",
#       class == "Novel Assembly" ~ "b",
#       class == "Novel Junction" ~ "c"
#     ),
#     class = factor(class, levels = levels(box_df$class)),
#     y = max(box_df$short_tpm, na.rm = TRUE) * 1.8
#   )

ggplot(box_df, aes(x = class, y = short_tpm, fill = class)) +
  geom_violin(
    alpha = 0.5,
    trim = TRUE,
    scale = "width"
  ) +
  geom_boxplot(
    width = 0.1,
    alpha = 0.6,
    outlier.shape = NA,
    colour = "black"
  ) +
  # geom_text(
  #   data = cld_df,
  #   aes(
  #     x = class,
  #     y = y,
  #     label = .group
  #   ),
  #   inherit.aes = FALSE,
  #   size = 5
  # ) +
  scale_y_continuous(
    trans = "log1p",
    breaks = breaks_shared,
    labels = scales::comma_format()
  ) +
  scale_fill_manual(values = my_palette) +
  scale_x_discrete(labels = c(
    "Shared" = "Shared",
    "Novel Assembly" = "Novel assembly",
    "Novel Junction" = "Novel junction(s)"
  )) +
  coord_cartesian(
    ylim = c(0, max(box_df$short_tpm, na.rm = TRUE) * 1.8)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    )
  ) +
  labs(
    x = "",
    y = "Short-read expression"
  )

# # Save
# ggsave(
#  filename = here::here("figures", 
#                        "13_sr_unique.png"),
#  plot = combo,
#  width = 9, 
#  height = 4, 
#  dpi = 1000
# )

################################################################################
### (5) Load and compare events discovered
################################################################################

# Load the IOE files (Change paths to your actual files)
base_path <- here::here("data",
                        "isoseq",
                        "15_suppa2")
modes <- c("short", "long")

# Function to load a single IOE and label it
load_mode_ioe <- function(mode) {
  # We can just pick one tissue folder (e.g., gonad) since the events are identical
  dir_path <- file.path(base_path, paste0(mode, "_gonad"), "events")
  files <- dir_ls(dir_path, glob = "*.ioe")
  
  map_df(files, function(f) {
    etype <- str_extract(basename(f), "(?<=events_)[A-Z0-9]+")
    
    read_tsv(f, show_col_types = FALSE) %>%
      mutate(mode = mode, event_type = etype) %>%
      dplyr::select(mode, event_type, gene_id, event_id)
  })
}

# Map over all combinations and combine into one master dataframe
all_ioe_collapsed <- map_df(modes, load_mode_ioe) %>%
  distinct(mode, event_type, gene_id, event_id)

# Filter out genes only present in one GTF so we can assess only shared 
# genes
shared_ioe_collapsed <- all_ioe_collapsed %>%
  filter(gene_id %in% shared_genes)

# Check
cat("Genes in unfiltered ioe files:", length(unique(all_ioe_collapsed$gene_id)), "\n")
cat("Genes in filtered ioe files:", length(unique(shared_ioe_collapsed$gene_id)), "\n")

### There is a discrepancy here where the genes in filtered = 8,553
### but the genes in shared_genes = 11,506

# Part of this is because single isoform genes can't find any events
# Identify the genes that have shared TPMs but NO events in the IOE
missing_from_ioe <- setdiff(shared_genes, shared_ioe_collapsed$gene_id) %>%
  unique()
# 2953 genes missing
# 2953 + 8553
# = 11506

# Check the isoform counts for these missing genes in your TPM data
missing_isoform_stats <- tpm_all_long %>%
  filter(gene %in% missing_from_ioe) %>%
  group_by(gene, data_type) %>%
  summarize(n_iso = n_distinct(transcript), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = data_type,
    values_from = n_iso,
    values_fill = 0
  ) %>%
  dplyr::rename(
    n_iso_SR = short,
    n_iso_LR = long
  ) %>%
  mutate(
    SR_multi = n_iso_SR > 1,
    LR_multi = n_iso_LR > 1
  )
# Summarize the findings
table(missing_isoform_stats$SR_multi, missing_isoform_stats$LR_multi)
### 2591 genes have only 1 isoform (in both datasets collapsed)
### 362 have more than one isoform - so why not events in those?

multi_no_ioe <- missing_isoform_stats %>%
     filter(SR_multi | LR_multi) %>%
     pull(gene)

### I checked manually in IGV and found the non-single isoform
### genes had differences in isoforms that were too compelx for
### SUPPA2 (ex. double exon skipping)

# Create comparison table
ioe_comparison <- shared_ioe_collapsed %>%
  group_by(event_id, event_type, gene_id) %>%
  summarize(
    in_short = any(mode == "short"),
    in_long = any(mode == "long"),
    .groups = "drop"
  ) %>%
  mutate(match_status = case_when(
    in_short & in_long   ~ "Shared",
    in_short & !in_long  ~ "Short",
    !in_short & in_long  ~ "Long"
  ))

# View the Results
table(ioe_comparison$event_type, ioe_comparison$match_status)
table(ioe_comparison$match_status)

### Make individual graphs for each event type
### Note: combined A3 and A5 categories

# Process the IOE comparisons from your pipeline
ioe_cleaned <- ioe_comparison %>%
  filter(!event_type %in% c("AF", "AL")) %>% # Remove AF/AL
  mutate(event_type = case_when(
    event_type %in% c("A3", "A5") ~ "A3/A5",
    TRUE ~ event_type
  )) %>%
  mutate(match_status = factor(
    match_status, 
    levels = c("Long", "Shared", "Short"),
    labels = c("Long", "Shared", "Short")
  ))

### Get # of events in each dataset for table 
table(ioe_cleaned$event_type, ioe_cleaned$in_short)
table(ioe_cleaned$event_type, ioe_cleaned$in_long)

table(ioe_cleaned$match_status)

# Calculate summary counts for plotting
plot_data_ioe <- ioe_cleaned %>%
  group_by(match_status, event_type) %>%
  summarize(count = n(), .groups = "drop")

# Define the structural order for the subplots
target_order <- c("SE", "MX", "A3/A5", "RI")
plot_data_ioe <- plot_data_ioe %>%
  mutate(event_type = factor(event_type, levels = target_order))

# Barplot function
make_individual_row <- function(data, target_event, show_x_axis = FALSE) {
  sub_df <- filter(data, event_type == target_event)
  
  # Calculate the total count for this event type to get the proportion
  total_event_count <- sum(sub_df$count, na.rm = TRUE)
  
  # Calculate a dynamic y-axis headroom limit to fit the numbers comfortably above the bars
  y_max <- max(sub_df$count, na.rm = TRUE) * 1.25  
  
  p <- ggplot(sub_df, aes(x = match_status, y = count, fill = event_type)) +
    # Increased width to 0.85 to make bars wider with less space between them
    geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.3) +
    # Add text labels (percentage only) above each bar
    geom_text(
      aes(label = paste0(round((count / total_event_count) * 100, 1), "%")),
      vjust = -0.5,       
      size = 3.5,         
      color = "#2c3e50"   # Not bolded
    ) +
    scale_fill_manual(values = c(
      "SE"    = "#94CCEC", 
      "MX"    = "#ECB4A0", 
      "A3/A5" = "#9DBCB5", 
      "RI"    = "#8676BF"
    )[target_event]) +
    # Set y-axis limits to prevent text clipping
    scale_y_continuous(labels = scales::comma, limits = c(0, y_max), expand = c(0, 0)) +
    labs(
      y = NULL,  
      x = NULL
    ) +
    theme_classic() +
    theme(
      axis.title.y = element_text(face = "bold", size = 11, angle = 90, vjust = 0.5, color = "#2c3e50"),
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey95", linewidth = 0.5),
      legend.position = "none"
    )
  
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    # Rotate the x-axis labels 45 degrees and align them cleanly
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  }
  
  return(p)
}

# Generate subplots vertically
row_plots <- list()
for (i in seq_along(target_order)) {
  current_event <- target_order[i]
  is_last <- (i == length(target_order))
  
  row_plots[[current_event]] = make_individual_row(
    data = plot_data_ioe, 
    target_event = current_event, 
    show_x_axis = is_last
  )
}

# Wrap only the right-side vertical panels together
final_vertical_stack <- wrap_plots(row_plots, ncol = 1)

# Render layout
final_vertical_stack

################################################################################
### (6) Parse SUPPA2 results
################################################################################

# I ran SUPPA2 in variable mode with fuzzy junctions 0
# I also decreased the # of replicates that have to have exp to return a result

# suppa2 results
sup_sum <- read.delim2(here::here("data",
                                  "isoseq",
                                  "15_suppa2",
                                  "events_summary_combined.tsv"))

# ensure numeric types (sometimes 'nan' forces the whole column to character)
sup_sum$p_value <- as.numeric(as.character(sup_sum$p_value))
sup_sum$dPSI <- as.numeric(as.character(sup_sum$dPSI))

# parse gene id and join with tpm_shared
sup_sum_parsed <- sup_sum %>%
  # regex breakdown:
  # ^([^;]+)    -> Group 1 (gene_id): Start to first semicolon
  # ;[^:]+:     -> Skip the semicolon and the event code (e.g., A3:)
  # ([^:]+)     -> Group 2 (chromosome): Everything until the next colon
  # :           -> The colon after chromosome
  # (.*)        -> Group 3 (coords): Everything (greedy) until...
  # :           -> ...the LAST colon in the string
  # ([^:]+)$    -> Group 4 (strand): Everything after that last colon until the end
  extract(
    Event_ID, 
    into = c("gene", "chromosome", "coords", "strand"), 
    regex = "^([^;]+);[^:]+:([^:]+):(.*):([^:]+)$",
    remove = TRUE
  ) %>%
  dplyr::rename(event_type = Event_Type) %>%
  filter(gene %in% shared_genes) %>% # ONLY LOOK AT RESULTS FOR SHARED GENES
  filter(!event_type %in% c("AF", "AL")) # REMOVE AF and AL

################################################################################
### (7) Plot all dpsi before coloured by sig
################################################################################

glimpse(sup_sum_parsed)

set_theme(theme_classic() + theme(strip.background = element_rect(fill = 'grey90', color = NA)))

# Reshape from long to wide format
sup_sum_wide <- sup_sum_parsed %>%
  pivot_wider(
    id_cols = c(tissue, gene, event_type, chromosome, coords, strand),
    names_from = mode,
    values_from = c(dPSI, p_value)
  )

# Explicitly remove events not tested in both platforms
sup_sum_wide <- sup_sum_wide %>%
  filter(!is.na(dPSI_long) & !is.na(dPSI_short))

# Grouping 'Long Only' and 'Short Only' significance
sup_sum_wide <- sup_sum_wide %>%
  mutate(sig_status = case_when(
    (p_value_long < 0.05 & abs(dPSI_long) > 0.2) & 
      (p_value_short < 0.05 & abs(dPSI_short) > 0.2) ~ "Both Significant",
    (p_value_long < 0.05 & abs(dPSI_long) > 0.2) | 
      (p_value_short < 0.05 & abs(dPSI_short) > 0.2) ~ "One Platform Only",
    TRUE ~ "Not Significant"
  ))

# Define colors for FILL and BORDER
# We use Shape 21 which supports both 'fill' and 'color' (outline)
fill_values <- c(
  "Both Significant"  = "black",
  "One Platform Only" = "grey40", 
  "Not Significant"   = "white"
)

border_values <- c(
  "Both Significant"  = "black",
  "One Platform Only" = "grey40", 
  "Not Significant"   = "grey85" # Light grey outline
)

# Create the Plot
max_val <- max(c(abs(sup_sum_wide$dPSI_long), abs(sup_sum_wide$dPSI_short)), na.rm = TRUE)

# Update Tissue Order
sup_sum_wide <- sup_sum_wide %>%
  mutate(tissue = factor(tissue, levels = c("tail", "gonad")))

# Calculate pearson correlation per tissue
sup_sum_wide %>%
  group_by(tissue) %>%
  summarise(
    test = list(cor.test(
      dPSI_long,
      dPSI_short,
      method = "pearson"
    ))
  ) %>%
  mutate(
    r = map_dbl(test, ~ .x$estimate),
    p = map_dbl(test, ~ .x$p.value)
  )

# Add pearson correlation to plot
# Calculate overall Pearson correlation
cor_label <- sup_sum_wide %>%
  summarise(
    r = cor(
      dPSI_long,
      dPSI_short,
      use = "complete.obs",
      method = "pearson"
    ),
    n = sum(complete.cases(dPSI_long, dPSI_short))
  ) %>%
  mutate(
    label = paste0(
      "r = ", round(r, 2),
      "\nn = ", scales::comma(n)
    ),
    x = 0.3 * max_val,
    y = -0.65 * max_val
  )

p4_combined <- ggplot(sup_sum_wide,
                      aes(x = dPSI_long, y = dPSI_short)) +
  # Quadrant colours
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#377EB8", alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E41A1C", alpha = 0.2) +
  # Correlation annotation
  geom_text(
    data = cor_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  # Plot points (ordering keeps significant points on top)
  geom_point(
    data = sup_sum_wide %>%
      arrange(factor(sig_status,
                     levels = c("Not Significant",
                                "One Platform Only",
                                "Both Significant"))),
    aes(fill = sig_status, color = sig_status),
    shape = 21,
    size = 1,
    stroke = 0.5,
    alpha = 0.8
  ) +
  scale_fill_manual(values = fill_values) +
  scale_color_manual(values = border_values) +
  
  coord_cartesian(
    xlim = c(-max_val, max_val),
    ylim = c(-max_val, max_val)
  ) +
  labs(
    x = "Long-read ΔPSI",
    y = "Short-read ΔPSI",
    fill = "Concordance",
    color = "Concordance"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

p4_combined

# Overall Pearson Correlation and p-value
overall_cor <- cor.test(
  sup_sum_wide$dPSI_long,
  sup_sum_wide$dPSI_short,
  method = "pearson"
)

# Print overall results
message("=== Overall Correlation ===")
message(paste("r:", round(overall_cor$estimate, 3)))
message(paste("n:", sum(complete.cases(sup_sum_wide$dPSI_long, sup_sum_wide$dPSI_short))))
message(paste("p-value:", format(overall_cor$p.value, scientific = FALSE)))

# 2. Per-Tissue Pearson Correlation and p-value
tissue_cor <- sup_sum_wide %>%
  group_by(tissue) %>%
  summarise(
    test = list(cor.test(dPSI_long, dPSI_short, method = "pearson")),
    n = n()
  ) %>%
  mutate(
    r = map_dbl(test, ~ .x$estimate),
    p = map_dbl(test, ~ .x$p.value),
    p_formatted = format(p, scientific = FALSE)
  ) %>%
  select(tissue, r, n, p, p_formatted)

message("\n=== Per-Tissue Correlation ===")
print(tissue_cor)


################################################################################
### (8) Compare significant events discovered
################################################################################

### Make individual bar plots

# Data wrangling

# Filter for Significant results (p < 0.05 & |dPSI| > 0.2)
significant_sup_sum <- sup_sum_parsed %>%
  filter(!is.na(p_value)) %>%            
  filter(p_value < 0.05) %>%
  filter(abs(dPSI) > 0.2)  

### Look at # genes in each
significant_sup_sum %>%
  group_by(tissue, mode) %>%
  summarise(
    n_genes = n_distinct(gene),
    .groups = "drop"
  )

# Filter and identify consistency across all tissues combined
plot_data_sig_match <- significant_sup_sum %>%
  group_by(gene, chromosome, event_type, coords) %>% 
  summarize(
    in_short = any(mode == "short"),
    in_long = any(mode == "long"),
    .groups = "drop"
  ) %>%
  mutate(match_status = case_when(
    in_short & in_long  ~ "Shared",
    in_short & !in_long ~ "Short",
    !in_short & in_long ~ "Long"
  )) %>%
  # NEW: Combine A3 and A5 into a unified event type prior to counting
  mutate(event_type = case_when(
    event_type %in% c("A3", "A5") ~ "A3/A5",
    TRUE ~ event_type
  )) %>%
  # Summarize for the plot
  group_by(match_status, event_type) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(
    match_status = factor(match_status, 
                          levels = c("Long", "Shared", "Short")),
    # Updated target structural order
    event_type = factor(event_type, levels = c("SE", "MX", "A3/A5", "RI"))
  )

# Update barplot function
make_individual_row_sig <- function(data, target_event, show_x_axis = FALSE) {
  sub_df <- filter(data, event_type == target_event)
  
  total_event_count <- sum(sub_df$count, na.rm = TRUE)
  y_max <- max(sub_df$count, na.rm = TRUE) * 1.25  
  
  p <- ggplot(sub_df, aes(x = match_status, y = count, fill = event_type)) +
    geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = paste0(round((count / total_event_count) * 100, 1), "%")),
      vjust = -0.5,       
      size = 3.5,         
      color = "#2c3e50"   
    ) +
    scale_fill_manual(values = c(
      "SE"    = "#94CCEC", 
      "MX"    = "#ECB4A0", 
      "A3/A5" = "#9DBCB5", 
      "RI"    = "#8676BF"
    )[target_event]) +
    scale_y_continuous(labels = scales::comma, limits = c(0, y_max), expand = c(0, 0)) +
    labs(
      y = NULL,  
      x = NULL
    ) +
    theme_classic() +
    theme(
      # FIX: Force the y axis title structure to be completely empty
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey95", linewidth = 0.5),
      legend.position = "none"
    )
  
  if (!show_x_axis) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  }
  
  return(p)
}

# Stack subplots
target_order <- c("SE", "MX", "A3/A5", "RI")
row_plots <- list()

for (i in seq_along(target_order)) {
  current_event <- target_order[i]
  is_last <- (i == length(target_order))
  
  # FIX: Updated to call *_sig variant function and matching data dataframe
  row_plots[[current_event]] = make_individual_row_sig(
    data = plot_data_sig_match, 
    target_event = current_event, 
    show_x_axis = is_last
  )
}

# Wrap panels together
final_vertical_stack_sig <- wrap_plots(row_plots, ncol = 1)

# Render layout window
final_vertical_stack_sig


### Find the # of genes that show sex-biased splicing in each tissue+mode
glimpse(significant_sup_sum)

table(significant_sup_sum$mode, significant_sup_sum$event_type)

gene_counts <- significant_sup_sum %>%
  group_by(tissue, mode) %>%
  summarise(
    n_unique_genes = n_distinct(gene),
    .groups = "drop"
  )

gene_counts

################################################################################
### (9) UpSet plots of events and sig events
################################################################################

glimpse(sup_sum_parsed)
glimpse(significant_sup_sum)

library(UpSetR)

# Set event order
event_order <- c("SE", "MX", "A3/A5", "RI")
event_order_upset <- rev(event_order)

### LONG - ALL EVENTS

# Filter and combine event types
long_all <- sup_sum_parsed %>%
  filter(mode == "long") %>%
  mutate(
    event_type = if_else(event_type %in% c("A3", "A5"),
                         "A3/A5",
                         event_type)
  ) %>%
  distinct(gene, event_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = event_type,
    values_from = value,
    values_fill = 0
  )

UpSetR::upset(
  as.data.frame(long_all[, event_order_upset]),
  sets = event_order_upset,
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = c(
    "#8676BF",  # RI
    "#9DBCB5",  # A3/A5
    "#ECB4A0",  # MX
    "#94CCEC"   # SE
  )
)

### SHORT - ALL EVENTS

# Filter and combine event types
short_all <- sup_sum_parsed %>%
  filter(mode == "short") %>%
  mutate(
    event_type = if_else(event_type %in% c("A3", "A5"),
                         "A3/A5",
                         event_type)
  ) %>%
  distinct(gene, event_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = event_type,
    values_from = value,
    values_fill = 0
  )

UpSetR::upset(
  as.data.frame(short_all[, event_order_upset]),
  sets = event_order_upset,
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = c(
    "#8676BF",  # RI
    "#9DBCB5",  # A3/A5
    "#ECB4A0",  # MX
    "#94CCEC"   # SE
  )
)

### LONG - ALL EVENTS

# Filter and combine event types
long_sig <- significant_sup_sum %>%
  filter(mode == "long") %>%
  mutate(
    event_type = if_else(event_type %in% c("A3", "A5"),
                         "A3/A5",
                         event_type)
  ) %>%
  distinct(gene, event_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = event_type,
    values_from = value,
    values_fill = 0
  )

UpSetR::upset(
  as.data.frame(long_sig[, event_order_upset]),
  sets = event_order_upset,
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = c(
    "#8676BF",  # RI
    "#9DBCB5",  # A3/A5
    "#ECB4A0",  # MX
    "#94CCEC"   # SE
  )
)

### SHORT - ALL EVENTS

# Filter and combine event types
short_sig <- significant_sup_sum %>%
  filter(mode == "short") %>%
  mutate(
    event_type = if_else(event_type %in% c("A3", "A5"),
                         "A3/A5",
                         event_type)
  ) %>%
  distinct(gene, event_type) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = event_type,
    values_from = value,
    values_fill = 0
  )

UpSetR::upset(
  as.data.frame(short_sig[, event_order_upset]),
  sets = event_order_upset,
  keep.order = TRUE,
  order.by = "freq",
  sets.bar.color = c(
    "#8676BF",  # RI
    "#9DBCB5",  # A3/A5
    "#ECB4A0",  # MX
    "#94CCEC"   # SE
  )
)

################################################################################
### (9) Stacked bar plot of sex bias in long-reads
################################################################################

glimpse(sup_sum_wide)

# Filter and categorize the data
plot_data <- sup_sum_wide %>%
  # Remove rows where dPSI_long is NaN or NA
  filter(!is.na(dPSI_long)) %>%
  # Filter for significant events only
  filter(sig_status != "Not Significant") %>%
  # Create a new column for Male vs Female bias
  mutate(bias = ifelse(dPSI_long > 0, "Female-biased", "Male-biased"))

# Create the faceted stacked bar plot
ggplot(plot_data, aes(x = bias, fill = event_type)) +
  geom_bar(position = "stack", color = "black", size = 0.2) +
  # Create separate plots for each tissue
  # scales = "free_y" allows the Y-axis to adjust to the count in each tissue
  facet_wrap(~tissue, scales = "free_y") + 
  scale_fill_brewer(palette = "Set3") + 
  labs(
    title = "Sex-Biased Splicing Events by Tissue",
    subtitle = "Significant events only (dPSI_long)",
    x = "Sex Bias Direction",
    y = "Number of Events",
    fill = "Event Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

### Redo with proportion

# Filter and categorize the data
plot_data <- sup_sum_wide %>%
  # Remove rows where dPSI_long is NaN or NA
  filter(!is.na(dPSI_long)) %>%
  # Filter for significant events only
  filter(sig_status != "Not Significant") %>%
  # Create a new column for Male vs Female bias
  mutate(bias = ifelse(dPSI_long > 0, "Female-biased", "Male-biased"))

# Create the proportional stacked bar plot
ggplot(plot_data, aes(x = bias, fill = event_type)) +
  # position = "fill" converts counts to proportions (0 to 1)
  geom_bar(position = "fill", color = "black", size = 0.2) +
  facet_wrap(~tissue) + 
  scale_fill_brewer(palette = "Set3") + 
  # Use scale_y_continuous to format the axis as percentages
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Sex-Biased Splicing Events",
    subtitle = "Significant events (dPSI_long) normalized by category",
    x = "Sex Bias Direction",
    y = "Percentage of Events",
    fill = "Event Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

library(dplyr)
library(ggplot2)

# 1. Filter and categorize the data using the corrected logic
# Added filter to keep only MX, SE, and RI event types
plot_data <- sup_sum_wide %>%
  filter(!is.na(dPSI_long)) %>%
  filter(sig_status != "Not Significant") %>%
  filter(event_type %in% c("MX", "SE", "RI")) %>% 
  mutate(bias = ifelse(dPSI_long > 0, "Female-biased", "Male-biased"))

# 2. Create the dodged bar plot
ggplot(plot_data, aes(x = event_type, fill = bias)) +
  # position = "dodge" places bars side-by-side
  geom_bar(position = "dodge", color = "black", linewidth = 0.2) +
  # Facet by tissue with free Y axes to see patterns in smaller tissues
  facet_wrap(~tissue, scales = "free_y") + 
  # Using distinct colors for sex-bias
  scale_fill_manual(values = c("Female-biased" = "#e41a1c", "Male-biased" = "#377eb8")) +
  labs(
    title = "Comparison of Sex-Biased Splicing Events by Tissue",
    subtitle = "Absolute counts of significant MX, SE, and RI events",
    x = "Event Type",
    y = "Number of Events",
    fill = "Bias Direction"
  ) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom"
  )

