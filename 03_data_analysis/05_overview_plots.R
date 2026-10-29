# Overview plots for first part of manuscript

# (1) Load
# (2) Isoform count histogram/stacked barplot
# (3) Isoform count stacked barplot with coding potential
# (4) Plot expression vs. potential and perplexity
# (5) Expression PCA
# (6) Clustered heatmap based on isoform expression
# (7) Ineffective vs. effective

# --- Data Manipulation & Wrangling ---
library(dplyr)
library(purrr)
library(reshape2)
library(stringr)
library(tidyr)
library(tidyverse) # Includes dplyr, ggplot2, tidyr, purrr, stringr

# --- Visualization & Formatting ---
library(circlize)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(ggplotify)
library(grid)
library(gt)
library(patchwork)
library(scales)
library(viridis)
library(gridExtra)

# --- Statistics ---
library(emmeans)
library(multcompView)
library(rstatix)
library(MASS)
library(broom)
library(ggpmisc)

################################################################################
### (1) Load
################################################################################

# Load mastertable
source("load_mt.R")

# Load sqanti output and orf status
sqanti_orf <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "16_edgeR_filterexp",
                                 "sqanti_w_orf.csv"))[,-1]

# Clean up
rm(mt_tissue, mt_samp)

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

# # Save
# ggsave(
#   filename = here::here("figures", "09_isoform_count_histogram.png"),
#   plot = iso_histogram,
#   width = 8,
#   height = 5,
#   dpi = 1000
# )

# Make histogram of relative abundance
mt_all %>%
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

# Plot
stacked_coding <- ggplot(plot_data2, 
                         aes(x = isoform_bin, 
                             y = proportion, 
                             fill = factor(iso_type, levels = level_order2))) +
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

# # Save
# ggsave(
#   filename = here::here("figures", "09_potential_histo_w_codingstatus.png"),
#   width = 8,
#   height = 5,
#   dpi = 1000
# )

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
  geom_point(color = "grey40", alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                 after_stat(rr.label), 
                                 sep = "~~~")),
               formula = y ~ x,
               parse = TRUE,
               label.x = "left",  # Fixed parameter name
               label.y = "top") +   # Fixed parameter name
  labs(
    x = "log2(TPM+1)",
    y = "Effective isoforms/Gene"
  ) +
  theme_minimal()


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

### COMBINE ALL TOGETHER

# Combine all together
combined_layout <- (
  iso_histogram + stacked_coding + plot_layout(guides = "collect")
) / 
  (p_potential + p_perplexity + p_tau + plot_layout(ncol = 3)) +
  plot_annotation(tag_levels = "A")  # use plot_annotation for overall labels

combined_layout

# # Save combined figure
# ggsave(
#   filename = here::here("figures", 
#                         "09_panel1.png"),
#   plot = combined_layout,
#   width = 12,   # adjust width and height as needed
#   height = 8,
#   dpi = 1000
# )

#################################################################################
### (5) Expression PCA
################################################################################

mt_samp <- read.csv(here::here("data",
                               "isoseq",
                               "11_mastertable",
                               "mastertable_sample.csv"))

### PCA based on gene expresseion
gene_mat_long <- mt_samp %>%
  group_by(sample, gene) %>%
  summarise(log2_geneTPM = mean(log2_geneTPM, na.rm = TRUE),
            .groups = "drop")

gene_mat_wide <- gene_mat_long %>%
  pivot_wider(
    names_from = gene,
    values_from = log2_geneTPM
  )

gene_mat <- gene_mat_wide %>%
  column_to_rownames("sample") %>%
  as.matrix()

gene_mat[is.na(gene_mat)] <- 0

gene_mat <- gene_mat[, apply(gene_mat, 2, var, na.rm = TRUE) > 0]

pca_res <- prcomp(gene_mat,
                  center = TRUE,
                  scale. = TRUE)

sample_info <- data.frame(
  sample = rownames(pca_res$x)
)

sample_info$sex <- ifelse(grepl("F", sample_info$sample),
                          "female", "male")

sample_info$tissue <- ifelse(grepl("G", sample_info$sample),
                             "gonad", "tail")

pca_df <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(sample_info, by = "sample")

# Percent variance explained
percent_var <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

# Create axis labels
pc1_lab <- paste0("PC1 (", round(percent_var[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(percent_var[2], 1), "%)")

p2 <- ggplot(pca_df, aes(PC1, PC2, color = sex, shape = tissue)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    title = "Gene Expression")


### PCA based on isoform expression
isoform_mat_long <- mt_samp %>%
  group_by(sample, transcript) %>%
  summarise(log2_isoformTPM = mean(log2_isoformTPM, na.rm = TRUE),
            .groups = "drop")

isoform_mat_wide <- isoform_mat_long %>%
  pivot_wider(
    names_from = transcript,
    values_from = log2_isoformTPM
  )

isoform_mat <- isoform_mat_wide %>%
  column_to_rownames("sample") %>%
  as.matrix()

isoform_mat[is.na(isoform_mat)] <- 0

isoform_mat <- isoform_mat[, apply(isoform_mat, 2, var, na.rm = TRUE) > 0]

pca_res <- prcomp(isoform_mat,
                  center = TRUE,
                  scale. = TRUE)

sample_info <- data.frame(
  sample = rownames(pca_res$x)
)

sample_info$sex <- ifelse(grepl("F", sample_info$sample),
                          "female", "male")

sample_info$tissue <- ifelse(grepl("G", sample_info$sample),
                             "gonad", "tail")

pca_df <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(sample_info, by = "sample")

# Percent variance explained
percent_var <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

# Create axis labels
pc1_lab <- paste0("PC1 (", round(percent_var[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(percent_var[2], 1), "%)")

p3 <- ggplot(pca_df, aes(PC1, PC2, color = sex, shape = tissue)) +
  geom_point(size = 2) +
  theme_classic() +
  labs(x = pc1_lab,
       y = pc2_lab,
       title = "Isoform Expression")
p2
p3

p23 <- p2/p3
p23

# Save
ggsave(
  filename = here::here("figures", "09_exp_pca.png"),
  width = 4,
  height = 6,
  dpi = 300
)

################################################################################
### (6) Clustered heatmap based on isoform expression
################################################################################

# Filter top variable isoforms
isoform_mat_t <- t(isoform_mat)
isoform_var <- apply(isoform_mat_t, 1, var)          # variance per isoform
top_isoforms <- names(sort(isoform_var, decreasing = TRUE))[1:1000]
isoform_mat_top <- isoform_mat_t[top_isoforms, , drop = FALSE]  # subset rows

#Simplify sample names (remove X and P)
colnames(isoform_mat_top) <- gsub("^X|P$", "", colnames(isoform_mat_top))

# Column annotations for sex + tissue
sample_info <- data.frame(
  sample = colnames(isoform_mat_top)
) %>%
  mutate(
    sex = factor(ifelse(grepl("F", sample), "Female", "Male")),
    tissue = factor(ifelse(grepl("G", sample), "Gonad", "Tail"))
  ) %>%
  column_to_rownames("sample")

# Annotation colors (user can customize tissue)
ann_colors <- list(
  sex = c(Female = "#E87D72", Male = "#55BCC2"),  # example red/blue
  tissue = c(Gonad = "black", Tail = "grey")  # you can change these
)

# Precompute column clustering
col_dist <- dist(t(isoform_mat_top), method = "euclidean")
col_hclust <- hclust(col_dist, method = "complete")

# Flip the dendrogram node manually
# This example flips the second branch (swap left/right)
# You can inspect col_hclust$labels and adjust order manually
col_hclust$order <- rev(col_hclust$order)  # simple full flip
# Or for selective flip, swap segments of col_hclust$order manually

# Draw heatmap with precomputed cluster
p4 <- pheatmap(
  isoform_mat_top,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = col_hclust,   # use precomputed and flipped dendrogram
  annotation_col = sample_info,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 10,
  angle_col = "45")

p4

# Convert it to a ggplot object
p4_gg <- as.ggplot(p4)

ggsave(
  filename = here::here("figures", "09_exp_isoformheatmap.png"),
  plot = p4_gg,
  width = 6.5,
  height = 8,
  dpi = 300
)

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
# ggsave(
#   filename = here::here("figures", "09_exp_alltogether.png"),
#   plot = together,
#   width = 10,  # Increased width to give the wider p4_gg room to breathe
#   height = 7,
#   dpi = 1000
# )

################################################################################
### (7) Ineffective vs. effective
################################################################################

mt_treat <- read.csv(here::here("data",
                               "isoseq",
                               "11_mastertable",
                               "mastertable_treatment.csv"))

# Load sqanti output and orf status
sqanti_orf <- read.csv(here::here("data",
                                  "01_build_longread_transcriptome",
                                  "16_edgeR_filterexp",
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
  theme_minimal(base_size = 14) +
    theme(
      # White backgrounds everywhere
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      # Facet strips: white with no border
      strip.background = element_rect(fill = "grey90", colour = NA),
      # Don't set strip text size
      strip.text = element_text(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

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

# # Save with a wider aspect ratio to accommodate the new hero panel
# ggsave(
#   here::here("figures", "09_RA_Ineffective_ManuscriptPlot.png"), 
#   final_mega_figure, 
#   width = 11, 
#   height = 9, 
#   dpi = 300
# )

