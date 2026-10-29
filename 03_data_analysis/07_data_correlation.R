# Plot mapping results from kallisto after correction in edgeR for long and short
# reads

# (1) Load data
# (2) Calculate correlation between read types + plot

library(dplyr)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(tidyr)
library(patchwork)
library(here)

################################################################################
#### (1) Load data
################################################################################

# Load long-read tpm
long_read <- read.csv(here::here("data",
                              "isoseq",
                              "11_mastertable",
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

# Put together and then calculate total gene expression
# (Keep 0s for calculating mean, as in mastertable script for long-reads)
short_read <- bind_rows(
  short_gonad %>% 
    dplyr::select(-X, -transcript) %>% 
    pivot_longer(-gene, names_to = "condition", values_to = "TPM"),
  
  short_tail %>% 
    dplyr::select(-X, -transcript) %>% 
    pivot_longer(-gene, names_to = "condition", values_to = "TPM")
) %>%
  group_by(gene) %>%
  summarize(
    mean_expression = mean(TPM, na.rm = TRUE)
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
message(paste("Pearson r:", round(pearson_r, 3), "| p-value:", format.pval(pearson_p, digits = 3)))
message(paste("Spearman rho:", round(spearman_rho, 3), "| p-value:", format.pval(spearman_p, digits = 3)))

# Plot using hexagonal density binning with trendlines
ggplot(combined_tpm, aes(x = log_long_tpm, y = log_short_tpm)) +
  geom_hex(bins = 50) +  
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 0.5) +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 1) +
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
