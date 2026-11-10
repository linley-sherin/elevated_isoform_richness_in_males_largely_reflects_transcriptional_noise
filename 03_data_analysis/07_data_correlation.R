## Step 07: Plot mapping results from kallisto after correction in edgeR for long and short
## reads
## Last modified by Linley Sherin 2026-11-03
#
# Runs locally.
# (1) Load data
# (2) Calculate correlation between read types + plot
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/07_data_correlation.R
#
# Reads:  data/03_data_analysis/03_mastertable/mastertable_all.csv
#         data/02_build_shortread_transcriptome/05_edgeR/
#           01_gonadtran_TMMTPMmatrix_short.csv, 02_tailtran_TMMTPMmatrix_short.csv
# Writes: figures/07_data_corr.png                                Fig 1C
#         data/03_data_analysis/07_data_correlation/plots/fig1C_corr.rds  Fig 1C panel for 14
#
# Overwrites its outputs every time it runs.

library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(here)

################################################################################
#### (1) Load data
################################################################################

# Load long-read tpm
long_read <- read.csv(here::here("data",
                              "03_data_analysis",
                              "03_mastertable",
                              "mastertable_all.csv")) %>%
  dplyr::select(gene, mean_geneTPM_all) %>%
  distinct() %>%
  dplyr::rename(long_tpm = mean_geneTPM_all)

# Load short-read tpm
short_gonad <- read.csv(here::here("data",
                                   "02_build_shortread_transcriptome",
                                   "05_edgeR",
                                   "01_gonadtran_TMMTPMmatrix_short.csv"))
short_tail <- read.csv(here::here("data",
                                  "02_build_shortread_transcriptome",
                                  "05_edgeR",
                                  "02_tailtran_TMMTPMmatrix_short.csv"))

# Number of samples (both tissues)
n_all <- ncol(short_gonad) + ncol(short_tail) - 6  # minus row name, transcript, gene in each

# Put together and then calculate total gene expression
# Gene TPM per sample = sum of its isoforms; mean over all samples, with 0 where
# the gene is not expressed in a tissue (same as mean_geneTPM_all for long reads)
short_read <- bind_rows(
  short_gonad %>%
    dplyr::select(-X, -transcript) %>%
    pivot_longer(-gene, names_to = "condition", values_to = "TPM"),

  short_tail %>%
    dplyr::select(-X, -transcript) %>%
    pivot_longer(-gene, names_to = "condition", values_to = "TPM")
) %>%
  group_by(gene, condition) %>%
  summarize(gene_tpm = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  group_by(gene) %>%
  summarize(
    mean_expression = sum(gene_tpm) / n_all
  ) %>%
  dplyr::rename(short_tpm = mean_expression)

rm(short_gonad, short_tail)

### Combine keeping only genes exp in both datasets
combined_tpm <- long_read %>%
  inner_join(short_read, by = "gene") %>%
  mutate(
    log_long_tpm = log2(long_tpm + 1),
    log_short_tpm = log2(short_tpm + 1)
  )

################################################################################
#### Calculate correlation between read types and plot
################################################################################

# Calculate correlation coefficients and p-value
pearson_test <- cor.test(combined_tpm$log_long_tpm,
                         combined_tpm$log_short_tpm,
                         method = "pearson")

# Calculate Spearman correlation and p-value
spearman_test <- cor.test(combined_tpm$log_long_tpm,
                          combined_tpm$log_short_tpm,
                          method = "spearman")

# Extract coefficients and p-values
pearson_r <- pearson_test$estimate
pearson_p <- pearson_test$p.value

spearman_rho <- spearman_test$estimate
spearman_p <- spearman_test$p.value

# Print results
message(paste("Genes in both datasets:", nrow(combined_tpm)))
message(paste("Pearson r:", round(pearson_r, 3), "| p-value:", format.pval(pearson_p, digits = 3)))
message(paste("Spearman rho:", round(spearman_rho, 3), "| p-value:", format.pval(spearman_p, digits = 3)))

# Plot using hexagonal density binning with trendlines
p_corr <- ggplot(combined_tpm, aes(x = log_long_tpm, y = log_short_tpm)) +
  geom_hex(bins = 50) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 1) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0("r = ", round(pearson_r, 3)),
    hjust = 1.1,
    vjust = -0.5,
    size = 4
  ) +
  scale_fill_viridis_c(
    option = "mako",
    trans = "log10",
    begin = 0.25,
    end = 0.75
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Long-read expression",
    y = "Short-read expression",
    fill = "Gene Count",
  )

p_corr

# Save
ggsave(here::here("figures", "07_data_corr.png"),
       plot = p_corr, width = 5, height = 4, dpi = 300)

plot_dir <- here::here("data", "03_data_analysis", "07_data_correlation", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(p_corr, file.path(plot_dir, "fig1C_corr.rds"))
