## Step 05: Overview plots for first part of manuscript
## Last modified by Linley Sherin 2026-11-03
#
# Runs locally.
# (1) Load
# (2) Isoform count histogram/stacked barplot
# (3) Isoform count stacked barplot with coding potential
# (4) Plot expression vs. potential and perplexity
# (5) Expression PCA
# (6) Clustered heatmap based on isoform expression
# (7) Ineffective vs. effective
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/05_overview_plots.R
#
# Reads:  data/03_data_analysis/03_mastertable/ (via lib/load_mt.R)
#         data/03_data_analysis/02_orf_status/sqanti_w_orf.csv
#         data/02_build_shortread_transcriptome/05_edgeR/0[12]_*_TMMTPMmatrix_short.csv
# Writes: data/03_data_analysis/05_overview_plots/plots/fig2[BCD]_*.rds   panels for 14_assemble_figures.R
#         figures/05_isoform_count_histogram.png      Fig 2B
#         figures/05_relabund_histogram.png           Fig 2C
#         figures/05_potential_histo_w_codingstatus.png   Fig 2D
#         figures/05_expression_vs_richness_diversity.png Supp Fig S5
#         figures/05_exp_pca.png, 05_exp_isoformheatmap.png, 05_exp_alltogether.png   Supp Fig S4
#         figures/05_exp_pca_short.png, 05_exp_isoformheatmap_short.png   short-read versions (reviewer)
#         data/03_data_analysis/05_overview_plots/plots/figS4_*.rds   panels for 14_assemble_figures.R
#         figures/05_RA_Ineffective_ManuscriptPlot.png    Supp Fig S6
#
# Overwrites its outputs every time it runs.

# --- Data Manipulation & Wrangling ---
library(dplyr)
library(tidyr)
library(tidyverse) # Includes dplyr, ggplot2, tidyr, purrr, stringr

# --- Visualization & Formatting ---
library(pheatmap)
library(ggplot2)
library(ggplotify)
library(patchwork)
library(scales)

# --- Statistics ---
library(ggpmisc)

################################################################################
### (1) Load
################################################################################

# Load mastertable
source(here::here("03_data_analysis", "lib", "load_mt.R"))

# Clean up
rm(mt_tissue)

# Panels are also saved as R objects so 14_assemble_figures.R can build the paper figures
plot_dir <- here::here("data", "03_data_analysis", "05_overview_plots", "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
### (2) Isoform count histogram + relative abundance histogram
################################################################################

# Compute isoform counts per gene
gene_iso_counts <- mt_all %>%
  group_by(gene) %>%
  summarise(n_isoforms = n_distinct(transcript), .groups = "drop") %>%
  mutate(
    isoform_bin = case_when(
      n_isoforms >= 5 ~ "5+",
      TRUE ~ as.character(n_isoforms)
    ),
    isoform_bin = factor(isoform_bin, levels = c("1","2","3","4","5+"))
  )

# Count number of genes per bin + compute %
hist_data <- gene_iso_counts %>%
  group_by(isoform_bin) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  mutate(
    percent = n_genes / sum(n_genes)
  )

# Create x-axis labels with counts
hist_labels <- hist_data %>%
  mutate(label = paste0(isoform_bin, "\n(n=", n_genes, ")"))

x_labels_hist <- setNames(hist_labels$label, hist_labels$isoform_bin)

# Plot (y = percent, labels on x-axis)
iso_histogram <- ggplot(hist_data, aes(x = isoform_bin, y = percent)) +
  geom_col(fill = "grey70", width = 0.9) +
#  scale_x_discrete(labels = x_labels_hist) +
  labs(
    x = "Isoform richness",
    y = "Genes"
  ) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),   # Increases X-axis label ("Isoform richness")
    axis.title.y = element_text(size = 14),   # Increases Y-axis label ("Isoforms")
    axis.text.x = element_text(size = 12, lineheight = 0.9), # Optional: increases X tick text
    axis.text.y = element_text(size = 12),    # Optional: increases Y tick text
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

iso_histogram

# Mean, median, and 95% CI for isoform richness per gene
mean_isoforms <- mean(gene_iso_counts$n_isoforms)
median_isoforms <- median(gene_iso_counts$n_isoforms)

se_isoforms <- sd(gene_iso_counts$n_isoforms) / sqrt(nrow(gene_iso_counts))
ci_isoforms <- qt(0.975, df = nrow(gene_iso_counts) - 1) * se_isoforms

cat("Mean isoforms per gene:", round(mean_isoforms, 2), "\n")
cat("Median isoforms per gene:", median_isoforms, "\n")
cat(
  "95% CI for mean:",
  round(mean_isoforms - ci_isoforms, 2),
  "–",
  round(mean_isoforms + ci_isoforms, 2),
  "\n"
)

# Calculate mean, median, and 95% CI for perplexity per gene

perplexity_data <- mt_all %>%
  dplyr::select(gene, pGene_all) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(!is.na(pGene_all))

# Mean and median
mean_perplexity <- mean(perplexity_data$pGene_all)
median_perplexity <- median(perplexity_data$pGene_all)

# 95% CI for the mean
se_perplexity <- sd(perplexity_data$pGene_all) /
  sqrt(nrow(perplexity_data))

ci_perplexity <- qt(
  0.975,
  df = nrow(perplexity_data) - 1
) * se_perplexity

# Print results
cat(
  "Mean perplexity per gene:",
  round(mean_perplexity, 2),
  "\n"
)

cat(
  "Median perplexity per gene:",
  round(median_perplexity, 2),
  "\n"
)

cat(
  "95% CI for mean:",
  round(mean_perplexity - ci_perplexity, 2),
  "–",
  round(mean_perplexity + ci_perplexity, 2),
  "\n"
)

# Save
ggsave(
  filename = here::here("figures", "05_isoform_count_histogram.png"),
  plot = iso_histogram,
  width = 8,
  height = 5,
  dpi = 1000
)
saveRDS(iso_histogram, file.path(plot_dir, "fig2B_isoform_count_histogram.rds"))

# Make histogram of relative abundance
ra_histogram <- mt_all %>%
  ggplot(aes(x = ra_all)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))), 
    fill = "grey70", 
    color = "white",
    binwidth = 0.2,     # 5 bins of 20% each
    boundary = 0,       # Aligns them perfectly at 0, 0.2, 0.4, etc.
    closed = "right"    # Safely captures 1 inside the 80%-100% bin
  ) +
  labs(
    x = "Relative abundance",
    y = "Isoforms"
  ) +
  scale_x_continuous(
    labels = scales::percent,
    breaks = c(0, 0.5, 1.0)
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(xlim = c(0, 1)) + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),   # Increases X-axis label ("Isoform richness")
    axis.title.y = element_text(size = 14),   # Increases Y-axis label ("Isoforms")
    axis.text.x = element_text(size = 12, lineheight = 0.9), # Optional: increases X tick text
    axis.text.y = element_text(size = 12),    # Optional: increases Y tick text
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ra_histogram

# Save
ggsave(
  filename = here::here("figures", "05_relabund_histogram.png"),
  plot = ra_histogram,
  width = 8,
  height = 5,
  dpi = 1000
)
saveRDS(ra_histogram, file.path(plot_dir, "fig2C_relabund_histogram.rds"))

# Count the number of rows in each relative abundance bin
bin_counts <- mt_all %>%
  mutate(
    ra_bin = cut(
      ra_all, 
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
      include.lowest = TRUE, 
      labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
    )
  ) %>%
  count(ra_bin, name = "isoform_count") %>%
  # Calculate proportion column
  mutate(isoform_proportion = (isoform_count / sum(isoform_count)*100))

print(bin_counts)

################################################################################
### (3) Isoform count stacked barplot with coding potential
################################################################################

# Compute isoform type for mt_all
mt_all2 <- mt_all %>%
  dplyr::select(gene, transcript, ptc_status, orf_identity) %>%
  group_by(gene) %>%
  arrange(orf_identity, transcript) %>%  # reproducible ordering
  mutate(
    orf_first = !duplicated(orf_identity)
  ) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding" ~ "Non-coding",
      ptc_status == "intron_retention" ~ "Retained intron",
      ptc_status == "nmd" ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  )

iso_counts <- sort(table(mt_all2$iso_type), decreasing = TRUE)
iso_prop <- prop.table(iso_counts)

cat("Counts for iso_type (sorted):\n")
cat(paste0(
  names(iso_counts), ": ",
  iso_counts, " (",
  sprintf("%.1f%%", iso_prop * 100),
  ")"
), sep = "\n")

# Compute isoform counts per gene with bins 1–5+
plot_data2 <- mt_all2 %>%
  group_by(gene) %>%
  mutate(isoform_count = n_distinct(transcript)) %>%
  ungroup() %>%
  mutate(
    isoform_bin = case_when(
      isoform_count >= 5 ~ "5+",
      TRUE ~ as.character(isoform_count)
    ),
    isoform_bin = factor(isoform_bin, 
                         levels = c("1","2","3","4","5+"))
  ) %>%
  group_by(isoform_bin, iso_type) %>%
  summarise(
    status_gene_sum = sum(1 / isoform_count), 
    .groups = "drop"
  ) %>%
  group_by(isoform_bin) %>%
  mutate(
    proportion = status_gene_sum / sum(status_gene_sum)
  ) %>%
  ungroup()

# Create x-axis labels with totals per bin
bin_totals2 <- plot_data2 %>%
  group_by(isoform_bin) %>%
  summarise(total = round(sum(status_gene_sum)), .groups = "drop") %>%
  mutate(axis_label = paste0(isoform_bin, "\n(n=", total, ")"))

x_labels2 <- setNames(bin_totals2$axis_label, bin_totals2$isoform_bin)

# Define colors and legend labels for iso_types
# iso_colors <- c(
#   "Non-coding"              = "#D4C2C6",
#   "Retained intron" = "#858893",
#   "Premature STOP"   = "#BB6C87",
#   "UTR variant"     = "#F8DCA7",
#   "Unique ORF"   = "#97C7C0" 
# )

iso_colors <- c(
  "Non-coding" = "#B89F9A",
  "Retained intron" = "#B5B8C1",
  "Premature STOP" = "#A17B6A",
  "UTR variant" = "#EBCE95",
  "Unique ORF" = "#81A19C" 
)

legend_labels2 <- c(
  "Non-coding"              = "Non-coding",
  "Premature STOP"   = "Premature STOP",
  "Retained intron" = "Retained intron",
  "Unique ORF"      = "Unique ORF",
  "UTR variant"     = "UTR variant"
)

# Define stacking order
level_order2 <- c(
  "Non-coding",
  "Premature STOP",
  "Retained intron",
  "UTR variant",
  "Unique ORF"
)

# Set stacking order in the data (so the saved plot does not need level_order2)
plot_data2$iso_type <- factor(plot_data2$iso_type, levels = level_order2)

# Plot
stacked_coding <- ggplot(plot_data2, 
                         aes(x = isoform_bin, 
                             y = proportion, 
                             fill = iso_type)) +
  geom_bar(stat = "identity", width = .9) +
  #scale_x_discrete(labels = x_labels2) +
  scale_y_continuous(labels = percent,
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = iso_colors,
    breaks = level_order2,
    labels = legend_labels2,
    guide = guide_legend(reverse = TRUE)  # Reverses legend key order
  ) +
  labs(
    x = "Isoform richness",
    y = "Isoforms",
    fill = "Protein-coding potential"     # Adds the legend title
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 14),   # Increases X-axis label ("Isoform richness")
    axis.title.y = element_text(size = 14),   # Increases Y-axis label ("Isoforms")
    axis.text.x = element_text(size = 12, lineheight = 0.9), # Optional: increases X tick text
    axis.text.y = element_text(size = 12),    # Optional: increases Y tick text
    panel.grid.major.x = element_blank()
  )

stacked_coding

# Save
ggsave(
  filename = here::here("figures", "05_potential_histo_w_codingstatus.png"),
  plot = stacked_coding,
  width = 8,
  height = 5,
  dpi = 1000
)
saveRDS(stacked_coding, file.path(plot_dir, "fig2D_potential_codingstatus.rds"))

### Get % genes that alternatively splice

# Get the total number of unique genes
total_genes <- n_distinct(mt_all$gene)

# Get the isoform counts per gene
gene_isoform_counts <- mt_all %>% 
  group_by(gene) %>% 
  summarise(n = n_distinct(transcript), .groups = "drop")

# Calculate percentages
p_2plus <- (sum(gene_isoform_counts$n >= 2) / total_genes) * 100
p_5plus <- (sum(gene_isoform_counts$n >= 5) / total_genes) * 100
p_10plus <- (sum(gene_isoform_counts$n >= 10) / total_genes) * 100

# Results
cat("% Genes with 2+ isoforms:", round(p_2plus, 2), "%\n")
cat("% Genes with 5+ isoforms:", round(p_5plus, 2), "%\n")
cat("% Genes with 10+ isoforms:", round(p_10plus, 2), "%\n")

# Calculate isoform count per gene
iso_stats <- mt_all %>%
  group_by(gene) %>%
  summarise(n_isoforms = n_distinct(transcript), .groups = "drop")

# Print Mean and Median
cat("Mean isoforms per gene: ", mean(iso_stats$n_isoforms), "\n")
cat("Median isoforms per gene:", median(iso_stats$n_isoforms), "\n")

# Get mean and median perplexity as well
perp <- mt_all %>%
  dplyr::select(gene, pGene_all) %>%
  distinct()

# Print Mean and Median
cat("Mean isoforms per gene: ", mean(perp$pGene_all), "\n")
cat("Median isoforms per gene:", median(perp$pGene_all), "\n")

################################################################################
### (4) Plot expression vs. potential and perplexity
################################################################################

glimpse(mt_all)

mt_all_g <- mt_all %>%
  dplyr::select(gene, potential_all, pGene_all, tau_all, log2_mean_geneTPM_all) %>%
  distinct() %>%
  mutate(
    tau_all = ifelse(potential_all == 1, 1, tau_all) # Set tau = 1 if potential = 1 (rather than NA)
  )

# (1) Potential vs expression (UNLOGGED)
p_potential <- ggplot(mt_all_g, 
                      aes(x = log2_mean_geneTPM_all, 
                          y = potential_all)) +
  geom_point(color = "grey40", alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  stat_poly_eq(aes(label = after_stat(rr.label)), 
               formula = y ~ x,
               parse = TRUE,
               label.x = "left",  
               label.y = "top") +  
  scale_x_continuous(
    # Tell ggplot exactly where to place the ticks in log2 space
    breaks = log2(c(1, 100, 10000)),
    labels = scales::trans_format(function(x) 2^x, scales::comma_format(accuracy = 1))
  ) +
  labs(
    x = "Gene expression",
    y = "Isoform richness"
  ) +
  theme_minimal()
p_potential

# Potential vs expression
model_potential <- lm(
  potential_all ~ log2_mean_geneTPM_all,
  data = mt_all_g
)

summary(model_potential)$r.squared
summary(model_potential)$coefficients[2, 4]

# (2) Perplexity vs expression (matched styling)
p_perplexity <- ggplot(mt_all_g,
                       aes(x = log2_mean_geneTPM_all, 
                           y = pGene_all)) +
  geom_point(color = "grey40", alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  stat_poly_eq(aes(label = after_stat(rr.label)), # Keeps only the R2 label
               formula = y ~ x,
               parse = TRUE,
               label.x = "left",  
               label.y = "top") +  
  scale_x_continuous(
    # Tell ggplot exactly where to place the ticks in log2 space
    breaks = log2(c(1, 100, 10000)),
    # Tell ggplot how to format those back into raw numbers with commas (no decimals)
    labels = scales::trans_format(function(x) 2^x, scales::comma_format(accuracy = 1))
  ) +
  labs(
    x = "Gene expression", # Updated to reflect unlogged scale
    y = "Isoform diversity"
  ) +
  theme_minimal()
p_perplexity

model_perplexity <- lm(
  pGene_all ~ log2_mean_geneTPM_all,
  data = mt_all_g
)

cat(
  "Perplexity vs expression:\n",
  "R² =", round(summary(model_perplexity)$r.squared, 3),
  "\np =", format.pval(summary(model_perplexity)$coefficients[2, 4], digits = 3),
  "\n"
)

# # (3) Tau vs expression
# p_tau <- ggplot(mt_all_g, 
#                 aes(x = log2_mean_geneTPM_all, y = tau_all)) +
#   geom_point(color = "grey40", alpha = 0.5, size = 1) +
#   geom_smooth(method = "lm", color = "black", se = TRUE) +
#   stat_poly_eq(aes(label = after_stat(rr.label)), # Keeps only the R2 label
#                formula = y ~ x,
#                parse = TRUE,
#                label.x = "right",  
#                label.y = "bottom") +  
#   scale_x_continuous(
#     # Tell ggplot exactly where to place the ticks in log2 space
#     breaks = log2(c(1, 100, 10000)), 
#     # Tell ggplot how to format those back into raw numbers with commas (no decimals)
#     labels = scales::trans_format(function(x) 2^x, scales::comma_format(accuracy = 1))
#   ) +
#   labs(
#     x = "Gene expression", # Updated to reflect unlogged scale
#     y = "Tau"
#   ) +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, 1))  # limit y-axis without dropping points
# p_tau

# Save richness and diversity vs expression together
p_expression <- p_potential + p_perplexity

ggsave(
  filename = here::here("figures", "05_expression_vs_richness_diversity.png"),
  plot = p_expression,
  width = 8,
  height = 4,
  dpi = 1000
)

#################################################################################
### (5) Expression PCA
################################################################################

# PCA and heatmap are made the same way for long and short reads
# (reviewer asked for both). Values are log2(TPM + 1); gene TPM is the sum
# of its isoforms; anything not expressed in a sample is 0.

### Expression matrices (samples x genes or isoforms)

# Long reads (sample-level mastertable)
long_gene_mat <- mt_samp %>%
  group_by(sample, gene) %>%
  summarise(value = mean(log2_geneTPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = value) %>%
  column_to_rownames("sample") %>%
  as.matrix()

long_isoform_mat <- mt_samp %>%
  group_by(sample, transcript) %>%
  summarise(value = mean(log2_isoformTPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = transcript, values_from = value) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Short reads (TMM-TPM matrices from the short-read pipeline)
short_tpm <- bind_rows(
  read.csv(here::here("data",
                      "02_build_shortread_transcriptome",
                      "05_edgeR",
                      "01_gonadtran_TMMTPMmatrix_short.csv"))[, -1],
  read.csv(here::here("data",
                      "02_build_shortread_transcriptome",
                      "05_edgeR",
                      "02_tailtran_TMMTPMmatrix_short.csv"))[, -1]
) %>%
  pivot_longer(cols = matches("^R[FM][GT][0-9]$"),  # sample columns (RFG1 etc.)
               names_to = "sample",
               values_to = "isoformTPM") %>%
  filter(!is.na(isoformTPM))  # each tissue only has its own samples

short_gene_mat <- short_tpm %>%
  group_by(sample, gene) %>%
  summarise(value = log2(sum(isoformTPM) + 1), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = value) %>%
  column_to_rownames("sample") %>%
  as.matrix()

short_isoform_mat <- short_tpm %>%
  mutate(value = log2(isoformTPM + 1)) %>%
  dplyr::select(sample, transcript, value) %>%
  pivot_wider(names_from = transcript, values_from = value) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Set missing values to 0 and drop features with no variance
clean_mat <- function(mat) {
  mat[is.na(mat)] <- 0
  mat[, apply(mat, 2, var) > 0]
}

long_gene_mat <- clean_mat(long_gene_mat)
long_isoform_mat <- clean_mat(long_isoform_mat)
short_gene_mat <- clean_mat(short_gene_mat)
short_isoform_mat <- clean_mat(short_isoform_mat)

### PCA

plot_pca <- function(mat, title) {
  pca_res <- prcomp(mat,
                    center = TRUE,
                    scale. = TRUE)

  # Percent variance explained
  percent_var <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

  pca_df <- as.data.frame(pca_res$x) %>%
    tibble::rownames_to_column("sample") %>%
    mutate(sex = ifelse(grepl("F", sample), "female", "male"),
           tissue = ifelse(grepl("G", sample), "gonad", "tail"))

  ggplot(pca_df, aes(PC1, PC2, color = sex, shape = tissue)) +
    geom_point(size = 2) +
    theme_classic() +
    labs(x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
         y = paste0("PC2 (", round(percent_var[2], 1), "%)"),
         title = title)
}

p2 <- plot_pca(long_gene_mat, "Gene Expression (long reads)")
p3 <- plot_pca(long_isoform_mat, "Isoform Expression (long reads)")
p2_short <- plot_pca(short_gene_mat, "Gene Expression (short reads)")
p3_short <- plot_pca(short_isoform_mat, "Isoform Expression (short reads)")

p2
p3
p2_short
p3_short

p23 <- p2/p3
p23_short <- p2_short/p3_short

# Save
ggsave(
  filename = here::here("figures", "05_exp_pca.png"),
  plot = p23,
  width = 4,
  height = 6,
  dpi = 300
)
ggsave(
  filename = here::here("figures", "05_exp_pca_short.png"),
  plot = p23_short,
  width = 4,
  height = 6,
  dpi = 300
)

################################################################################
### (6) Clustered heatmap based on isoform expression
################################################################################

# Annotation colors (user can customize tissue)
ann_colors <- list(
  sex = c(Female = "#E87D72", Male = "#55BCC2"),  # example red/blue
  tissue = c(Gonad = "black", Tail = "grey")  # you can change these
)

plot_heatmap <- function(mat) {
  # Filter top variable isoforms
  mat_t <- t(mat)
  isoform_var <- apply(mat_t, 1, var)          # variance per isoform
  top_isoforms <- names(sort(isoform_var, decreasing = TRUE))[1:1000]
  mat_top <- mat_t[top_isoforms, , drop = FALSE]  # subset rows

  # Column annotations for sex + tissue
  sample_info <- data.frame(
    sample = colnames(mat_top)
  ) %>%
    mutate(
      sex = factor(ifelse(grepl("F", sample), "Female", "Male")),
      tissue = factor(ifelse(grepl("G", sample), "Gonad", "Tail"))
    ) %>%
    column_to_rownames("sample")

  # Precompute column clustering
  col_dist <- dist(t(mat_top), method = "euclidean")
  col_hclust <- hclust(col_dist, method = "complete")

  # Flip the dendrogram node manually
  # This example flips the second branch (swap left/right)
  # You can inspect col_hclust$labels and adjust order manually
  col_hclust$order <- rev(col_hclust$order)  # simple full flip
  # Or for selective flip, swap segments of col_hclust$order manually

  # Draw heatmap with precomputed cluster
  hm <- pheatmap(
    mat_top,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = col_hclust,   # use precomputed and flipped dendrogram
    annotation_col = sample_info,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,   # sample names too crowded to read in the figure
    fontsize = 10,
    angle_col = "45",
    silent = TRUE)

  # Convert it to a ggplot object
  as.ggplot(hm)
}

p4_gg <- plot_heatmap(long_isoform_mat)
p4_gg_short <- plot_heatmap(short_isoform_mat)

p4_gg
p4_gg_short

ggsave(
  filename = here::here("figures", "05_exp_isoformheatmap.png"),
  plot = p4_gg,
  width = 6.5,
  height = 8,
  dpi = 300
)
ggsave(
  filename = here::here("figures", "05_exp_isoformheatmap_short.png"),
  plot = p4_gg_short,
  width = 6.5,
  height = 8,
  dpi = 300
)

# Panels for Supp Fig S4 (built in 14_assemble_figures.R)
saveRDS(p2, file.path(plot_dir, "figS4_pca_gene_long.rds"))
saveRDS(p3, file.path(plot_dir, "figS4_pca_isoform_long.rds"))
saveRDS(p4_gg, file.path(plot_dir, "figS4_heatmap_long.rds"))
saveRDS(p2_short, file.path(plot_dir, "figS4_pca_gene_short.rds"))
saveRDS(p3_short, file.path(plot_dir, "figS4_pca_isoform_short.rds"))
saveRDS(p4_gg_short, file.path(plot_dir, "figS4_heatmap_short.rds"))

### PATCHWORK

# Create the left-hand column and collect its legends specifically
left_col <- (p2 / p3) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')

# Combine with p4_gg and set the relative widths
together <- (left_col | p4_gg) + 
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = 'A')

together

# Save
ggsave(
  filename = here::here("figures", "05_exp_alltogether.png"),
  plot = together,
  width = 10,  # Increased width to give the wider p4_gg room to breathe
  height = 7,
  dpi = 1000
)

################################################################################
### (7) Ineffective vs. effective
################################################################################

# Load sqanti output and orf status
sqanti_orf <- read.csv(here::here("data",
                                  "03_data_analysis",
                                  "02_orf_status",
                                  "sqanti_w_orf.csv"))[,-1]

glimpse(mt_treat)

# Global Preprocessing
mt_treat2 <- mt_treat %>%
  mutate(
    # Tail first, then other tissues (add others to the vector as needed)
    tissue = factor(
      tissue,
      levels = c("tail", "gonad"),
      labels = c("Tail", "Gonad")
    ),
    # Female first
    sex = factor(
      sex, 
      levels = c("female", "male"),
      labels = c("Female", "Male")),
    # Effective first
    effGene_treat = factor(effGene_treat, 
                           levels = c(TRUE, FALSE), 
                           labels = c("Effective", "Ineffective")))

# Consistent Styling
theme_analysis <- function() {
  theme_classic() +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text = element_text(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

sex_colors <- c("Female" = "#E87D72", "Male" = "#55BCC2")

# P1: Log2 Expression
p1 <- ggplot(mt_treat2, aes(x = effGene_treat, y = log2_mean_isoformTPM_treat, fill = sex)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  facet_grid(sex ~ tissue) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Expression level", x = NULL, y = "log2(mean TPM + 1)") +
  theme_analysis() +
  theme(legend.position = "none")

# P2: Variability
p2 <- ggplot(mt_treat2, aes(x = effGene_treat, y = cv_isoformTPM_treat, fill = sex)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
  facet_grid(sex ~ tissue) +
  scale_fill_manual(values = sex_colors) +
  labs(title = "Expression Variability", x = NULL, y = "Coefficient of Variation") +
  theme_analysis() +
  theme(legend.position = "none")

# Update P1 & P2 (Sex Legend)
p1 <- p1 + scale_fill_manual(values = sex_colors, name = "Sex", labels = c("Female", "Male"))
p2 <- p2 + scale_fill_manual(values = sex_colors, name = "Sex", labels = c("Female", "Male"))

# P3: PTC Status Composition
p3_data <- mt_treat2 %>%
  group_by(sex, tissue, effGene_treat, ptc_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sex, tissue, effGene_treat) %>%
  mutate(
    proportion = count / sum(count),
    ptc_status = recode(ptc_status, "coding"="Coding", "non-coding"="Non-coding", 
                        "nmd"="NMD", "intron_retention"="RI")
  )

# Put Coding LAST so it sits at the BOTTOM of the stack
ptc_levels <- c("Non-coding", "RI", "NMD", "Coding")
p3_data$ptc_status <- factor(p3_data$ptc_status, levels = ptc_levels)

# Since we want Turquoise (the 1st Brewer color) for Coding (the 4th level),
# we define a manual named vector to ensure the color stays linked to the name.
ptc_palette <- c(
  "Coding"      = "#66C2A5", # Turquoise (Set2 Color 1)
  "NMD"         = "#FC8D62", # Orange (Set2 Color 2)
  "RI"          = "#8DA0CB", # Blue-ish (Set2 Color 3)
  "Non-coding"  = "#E78AC3"  # Pink (Set2 Color 4)
)

p3 <- ggplot(p3_data, aes(x = effGene_treat, y = proportion, fill = ptc_status)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(sex ~ tissue) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Coding Status", x = NULL, y = "Proportion", fill = "PTC Status") +
  theme_analysis()

p3 <- p3 + scale_fill_manual(values = ptc_palette)

# P4: Annotation (Novel vs. Annotated)

# Label genes and isoforms based on whether they are novel or not
sqanti_orf <- sqanti_orf %>%
  mutate(
    annotated = ifelse(grepl("novel", 
                             associated_gene, 
                             ignore.case = TRUE),
                       "novel", 
                       "annotated"),
    isoform_category = ifelse(grepl("novel", associated_transcript, ignore.case = TRUE),
                              "novel", "annotated"))

# Combine with mt_treat
mt_treat2 <- mt_treat2 %>%
  left_join(dplyr::select(sqanti_orf, c(gene, transcript, annotated, isoform_category)))

p4_data <- mt_treat2 %>%
  group_by(sex, tissue, effGene_treat, isoform_category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sex, tissue, effGene_treat) %>%
  mutate(proportion = count / sum(count))

p4 <- ggplot(p4_data, aes(x = effGene_treat, y = proportion, fill = isoform_category)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(sex ~ tissue) +
  scale_fill_manual(values = c("annotated" = "grey70", "novel" = "#f4a1a1")) +
  labs(title = "Annotated vs. Novel", x = NULL, y = "Proportion", fill = "Category") +
  theme_analysis()

# Update P3 (PTC/Coding Status Legend)
# Note: We use ptc_palette which we defined earlier
p3 <- p3 + scale_fill_manual(values = ptc_palette, name = "Coding Status")


# Update P4 (Category Legend)
p4 <- p4 + scale_fill_manual(
  values = c("annotated" = "grey70", "novel" = "#f4a1a1"), 
  name = "Category",
  labels = c("Annotated", "Novel")
)

# Patchwork layout
# A is the hero panel (p5) on the left.
# B, C, D, E are the smaller panels on the right.
layout <- "
AABB
CCDD
"

# We list p5 first so it gets Tag A
final_mega_figure <- p1 + p2 + p3 + p4 + 
  plot_layout(
    design = layout, 
    guides = 'collect'
  ) + 
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = 'right')

# Display
final_mega_figure

# Save with a wider aspect ratio to accommodate the new hero panel
ggsave(
  here::here("figures", "05_RA_Ineffective_ManuscriptPlot.png"),
  final_mega_figure,
  width = 11,
  height = 9,
  dpi = 300
)

