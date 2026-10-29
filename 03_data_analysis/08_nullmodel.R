# Investigate data qualities

# (1) Load
# (2) Intron number distribution
# (3) Expression distribution
# (4) Compare input distributions across condition
# (5) Save for Wouter
# (6) Wouter's code - null model
# (7) Plot scatterplots
# (8) Plot quantiles

library(here)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(furrr)
library(GenomicRanges)
library(tibble)
library(scales)
library(tidyverse)
library(ggpubr)

################################################################################
### (1) Load
################################################################################

# Load master table for all samples
mt_treat <- read.csv(here::here("data",
                                "isoseq",
                                "11_mastertable",
                                "mastertable_treatment.csv"))

# Load metadata
meta <- read.csv(here::here("data",
                            "isoseq", 
                            "metadata.csv"))

# Load long-read transcriptome
ann_lg <- rtracklayer::import(here::here("data",
                                         "GTF",
                                         "04_filtered_by_exp_long.gtf"))

glimpse(ann_lg)

################################################################################
### (2) Intron number distribution
################################################################################

### Longest transcript approach

# Look only at exons
exons <- ann_lg[ann_lg$type == "exon"]

# Count exons per transcript
tx_exon <- as.data.frame(exons) |>
  group_by(gene_id, transcript_id) |>
  summarise(n_exons = n(), .groups = "drop")

# Pick longest transcript per gene
ref_tx <- tx_exon |>
  group_by(gene_id) |>
  slice_max(n_exons, n = 1, with_ties = FALSE) |>
  ungroup()

# Introns = exons - 1
ref_tx <- ref_tx |>
  mutate(n_introns = pmax(n_exons - 1, 0))

# Make sure it's only genes that made it through filtering
ref_tx <- ref_tx %>%
  filter(gene_id %in% mt_treat$gene)

# Histogram
ggplot(ref_tx, aes(x = n_introns)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(
    x = "Number of introns per gene",
    y = "Number of genes",
    title = "Distribution of intron counts"
  ) +
  theme_minimal()

################################################################################
### (3) Expression distribution
################################################################################

glimpse(mt_treat)

# We want to collapse by gene
exp <- mt_treat %>%
  group_by(tissue, sex) %>%
  dplyr::select(gene, mean_geneTPM_treat, log2_mean_geneTPM_treat, potential_treat, pGene_treat, tau_treat) %>%
  distinct() 

# Join with ref_tx
gene_feat <- exp %>%
  inner_join(ref_tx, by = c("gene" = "gene_id")) %>%
  dplyr::select(-transcript_id)

################################################################################
### (4) Compare input distributions across condition
################################################################################

gene_feat$tissue <- factor(
  gene_feat$tissue,
  levels = c("tail", "gonad"))

means <- gene_feat %>%
  group_by(sex, tissue) %>%
  summarise(mean_expr = mean(log2_mean_geneTPM_treat),
            mean_intron = mean(n_introns),
            .sexs = "drop")

# Expression distribution
exp_plot <- ggplot(
  gene_feat,
  aes(x = log2_mean_geneTPM_treat,
      fill = sex,
      color = sex)
) +
  geom_density(alpha = 0.2) +
  geom_vline(
    data = means,
    aes(xintercept = mean_expr,
        color = sex),
    linewidth = 0.5,
    linetype = 2
  ) +
  facet_wrap(~ tissue) +
  theme_bw()

# Intron distribution
intron_plot <- ggplot(
  gene_feat,
  aes(x = n_introns,
      fill = sex,
      color = sex)
) +
  geom_density(alpha = 0.2) +
  geom_vline(
    data = means,
    aes(xintercept = mean_intron,
        color = sex),
    linewidth = 0.5,
    linetype = 2
  ) +
  facet_wrap(~ tissue) +
  theme_bw()

# Generate correlation labels
cor_labels <- gene_feat %>%
  group_by(sex, tissue) %>%
  summarise(
    r2 = cor(n_introns, log2_mean_geneTPM_treat, use = "complete.obs")^2,
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(r2 < 0.01,
                   "R² < 0.01",
                   paste0("R² = ", sprintf("%.3f", r2)))
  )

# Plot against each other
versus_plot <- ggplot(gene_feat,
                      aes(x = n_introns,
                          y = log2_mean_geneTPM_treat)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  stat_regline_equation(
    aes(label = after_stat(eq.label)),
    size = 4
  ) +
  geom_text(
    data = cor_labels,
    aes(label = label),
    x = -Inf, y = Inf,
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_grid(sex ~ tissue) +
  theme_bw() +
  theme(
    text = element_text(size = 12)
  )

combined_plot <- ((exp_plot / intron_plot) | versus_plot) +
  plot_annotation(tag_levels = "A")

combined_plot

################################################################################
### (5) Save for Wouter
################################################################################

# Save
write.csv(gene_feat, here::here("data", 
                                "isoseq",
                                "13_nullmodel_wouter",
                                "gene_features_wouter.csv"),
          row.names = FALSE)

################################################################################
### (6) Wouter's code - null model
################################################################################

set_theme(theme_classic() + theme(strip.background = element_rect(fill = 'grey90', color = NA)))

#gene_feat <- data.table::fread('gene_features_wouter.csv') |> 
#  mutate(potentialisoforms = readr::parse_number(potentialisoforms))
gonad <- read.csv(here::here("data",
                              "01_build_longread_transcriptome",
                              "16_edgeR_filterexp",
                              "01_gonadtran_countmatrix.csv"))[,-1]
tail <- read.csv(here::here("data",
                              "01_build_longread_transcriptome",
                              "16_edgeR_filterexp",
                              "02_tailtran_countmatrix.csv"))[,-1]

# gonad: sum F columns and M columns
gonad_depths <- gonad %>%
  dplyr::select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>%
  mutate(
    sex = if_else(str_detect(sample, "M"), "male", "female"),
    tissue = "gonad"
  ) %>%
  group_by(tissue, sex) %>%
  summarise(total_reads = sum(count, na.rm = TRUE), .groups = "drop")

# tail: sum F columns and M columns
tail_depths <- tail %>%
  dplyr::select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "count") %>%
  mutate(
    sex = if_else(str_detect(sample, "M"), "male", "female"),
    tissue = "tail"
  ) %>%
  group_by(tissue, sex) %>%
  summarise(total_reads = sum(count, na.rm = TRUE), .groups = "drop")

all_depths <- bind_rows(gonad_depths, tail_depths)

# Join depths to your existing features
gene_feat_grouped <- gene_feat %>%
  inner_join(all_depths, by = c("tissue", "sex"))

# Core simulation functions
perplexity <- function(counts) { # counts is reads basically
  p <- counts / sum(counts)
  2 ^ (-sum(p * log2(p), na.rm = TRUE))
}

tau_specificity <- function(counts) {
  if (length(counts) <= 1) return(0)
  if (max(counts) == 0) return(0)
  x <- counts / max(counts)
  sum(1 - x) / (length(x) - 1)
}

simulate_gene <- function(n, n_introns, M) {
  if (n_introns == 0) { # If there is only 1 isoform
    return(data.frame(gene_potential = NA, gene_perplexity = NA, gene_tau = NA
    ))
  }
  splices <- matrix(rbinom(n_introns * n, 1, 1 - M), nrow = n_introns)
  isoforms <- apply(splices, 2, paste0, collapse = "")
  isoforms_tab <- table(isoforms)
  data.frame(
    gene_potential = length(isoforms_tab),
    gene_perplexity = perplexity(isoforms_tab),
    gene_tau = tau_specificity(isoforms_tab)
  )
}

# Wrapper to simulate a whole tissue-sex block
produce_one_dataset <- function(df, M) {
  # Calculate 'n' based on the specific group's read depth
  group_depth <- unique(df$total_reads)
  df_sim <- df %>%
    mutate(n = ceiling(mean_geneTPM_treat / sum(mean_geneTPM_treat) * group_depth)) %>%
    filter(n > 0)
  out <- future_map2_dfr(df_sim$n, df_sim$n_introns, simulate_gene, M = M, 
                         .options = furrr_options(seed = TRUE))
  bind_cols(dplyr::select(df_sim, gene, tissue, sex, mean_geneTPM_treat), out)
}

plan(multisession, workers = 24)

# run simulations
sim_data <- gene_feat_grouped %>%
  ungroup() %>%   # remove any existing grouping
  group_split(tissue, sex) %>%
  map_dfr(function(slice) {
    map_dfr(setNames(c(0, 0.007, 0.01), c(0, 0.007, 0.01)), 
            ~produce_one_dataset(slice, M = .x), .id = "error_rate")
  }) %>%
  mutate(type = "simulated")

# prep observed data
obs_data <- gene_feat_grouped %>%
  dplyr::select(
    gene, tissue, sex, mean_geneTPM_treat, 
    gene_perplexity = pGene_treat, 
    gene_potential = potential_treat,
    gene_tau = tau_treat  # just take it directly
  ) %>%
  mutate(
    type = "observed",
    error_rate = "observed"
  )

# combine and create a Group ID for coloring
full_results <- bind_rows(sim_data, obs_data) %>%
  mutate(
    group_id = paste(tissue, sex, sep = "_"),
    # Convert error_rate to factor to control facet order
    error_rate = factor(error_rate, levels = c("observed", "0", "0.007", "0.01"))
  )
plan(sequential)

# Sub tau = 1 for any gene where potential = 1
# (otherwise returns NA)
full_results <- full_results %>%
  mutate(
    gene_tau = if_else(gene_potential == 1, 1, gene_tau)
  )

################################################################################
### (7) Plot scatterplots
################################################################################

# Data prep
plot_data_final <- full_results %>%
  filter(error_rate %in% c("observed", "0.007")) %>%
  mutate(
    log2_tpm = log2(mean_geneTPM_treat + 1),
    group_id = factor(
      paste(sex, tissue, sep = "_"),
      levels = c("female_tail", "male_tail", "female_gonad", "male_gonad"),
      labels = c("Female Tail", "Male Tail", "Female Gonad", "Male Gonad")
    )
  )

# Plot function (updated to use lm)
plot_overlapping_comparison <- function(data, y_var, y_label,
                                        show_xlabel = FALSE,
                                        show_facets = FALSE) {
  ggplot(data, aes(x = log2_tpm, y = !!sym(y_var),
                   color = type, fill = type)) +
    geom_point(alpha = 0.1, size = 1.0, stroke = 0, show.legend = FALSE) +
    stat_density_2d(
      aes(alpha = after_stat(level)),
      geom = "polygon",
      bins = 6,
      adjust = 2,
      show.legend = FALSE
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1
    ) +
    stat_poly_eq(
      data = ~ subset(., type == "observed"),
      aes(
        label = after_stat(
          paste0(
            rr.label,
            "*\", \"*",
            ifelse(
              p.value < 0.001,
              "italic(p) < 0.001",
              paste0("italic(p) == ", format.pval(p.value, digits = 3, eps = 0))
            )
          )
        )
      ),
      formula = y ~ x,
      parse = TRUE,
      label.x = "left",
      label.y = "top",
      size = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~group_id, nrow = 1) +
    scale_fill_manual(values = c(
      "observed" = "#895293",
      "simulated" = "#BFB29D"
    )) +
    scale_color_manual(values = c(
      "observed" = "#895293",
      "simulated" = "#BFB29D"
    )) +
    scale_alpha_continuous(range = c(0.1, 0.4)) +
    labs(
      x = if (show_xlabel) "Gene Expression (log2(TPM + 1))" else NULL,
      y = y_label
    ) +
    theme_classic() +
    theme(
      strip.background = if (show_facets)
        element_rect(fill = "grey90", color = NA)
      else element_blank(),
      strip.text = if (show_facets)
        element_text(face = "plain")
      else element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.text.x = if (show_xlabel)
        element_text()
      else element_blank(),
      axis.ticks.x = if (show_xlabel)
        element_line()
      else element_blank()
    )
  }

# Plots with equations

p1 <- plot_overlapping_comparison(
  plot_data_final,
  "gene_potential",
  "Isoform\nrichness",
  show_xlabel = FALSE,
  show_facets = TRUE
) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000),
    labels = label_number()
  )

p2 <- plot_overlapping_comparison(
  plot_data_final,
  "gene_perplexity",
  "Isoform\ndiversity",
  show_xlabel = TRUE,
  show_facets = FALSE
) +
  coord_cartesian(ylim = c(1, 10))

# p3 <- plot_overlapping_comparison(
#   plot_data_final,
#   "gene_tau",
#   "Tau (Isoform Specificity)",
#   show_xlabel = TRUE,
#   show_facets = FALSE
# ) +
#   coord_cartesian(ylim = c(0.6, 1.0))

# Combine
final_plot <- (p1 / p2) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(legend.position = "right")

final_plot

# Save
ggsave(
  here::here("figures", "12_null_comparison.png"),
  plot = final_plot,
  width = 10,
  height = 8,
  dpi = 600
)
