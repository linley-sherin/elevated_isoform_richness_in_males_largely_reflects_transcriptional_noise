# Plots and stats for sex isoform comparisons

# (1) Load

### ISOFORM RICHNESS
# (2) Calculate effect size of sex for isoform richness
# (3) Leave one out analysis richness
# (4) Box plot for isoform richness (replicate means)
# (5) Orf barplots and forest plots (isoform richness)
# (6) Make log2FC for potential vs. bias (isoform richness) GONAD ONLY
# (7) Box plot for isoform richness vs bias GONAD ONLY
# (8) Isoform composition female biased vs. male biased GONAD ONLY

### ISOFORM DIVERSITY
# (9) Calculate effect size of sex for perplexity (isoform diversity)
# (10) Leave one out analysis diversity
# (11) Box plot for isoform diversity (replicate means)
# (12) Orf barplots and forest plots (isoform diversity)
# (13) Make log2FC for perplexity vs. bias (isoform diversity) GONAD ONLY
# (14) Box plot for isoform diversity vs bias GONAD ONLY
# (15) Eff isoform composition female biased vs. male biased GONAD ONLY

### BOTH
# (16) Bias orf forest plot (richness and diversity)

# --- Data Manipulation & Wrangling ---
library(purrr)
library(reshape2)
library(stringr)
library(fgsea)

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
library(gghalves)
library(cowplot)
library(ggnewscale)

# --- Statistics ---
library(emmeans)
library(multcompView)
library(rstatix)
library(MASS)
library(broom)

################################################################################
### (1) Load
################################################################################

# Load edgeR results
source("load_dge.R")

# Load master table
source("load_mt.R")

# Join with edgeR results and orf_status column of sqanti_w_orf
mt_treat <- mt_treat %>%
 left_join(combined_dge, by = c("gene", "tissue"))

# Clean up
rm(dge_gonad,
   dge_tail,
   dte_gonad,
   dte_tail,
   dtug_gonad,
   dtug_tail,
   dtut_gonad,
   dtut_tail,
   combined_dte,
   combined_dtug,
   combined_dtut,
   mt_tissue,
   mt_all)
 
################################################################################
### (2) Calculate effect size of sex for potential (isoform richness)
################################################################################

# Look at sample level
glimpse(mt_samp)
# Remember: All sex-specific genes already removed

### PUT BACK IN 0s FOR REPLICATES WITHOUT EXP
# Collapse to gene level and pivot wide, sub in 0 for NA
mt_ready <- mt_samp %>%
  dplyr::select(gene, sample, potential_samp) %>%
  distinct() %>%
  pivot_wider(
    names_from = sample, 
    values_from = potential_samp,
    values_fill = 0 
  )

# Split tissues safely with na.rm = TRUE
mt_gonad <- mt_ready %>%
  dplyr::select(gene, ends_with("GP")) %>%
  filter(rowSums(dplyr::select(., ends_with("GP")), na.rm = TRUE) > 0)

mt_tail <- mt_ready %>%
  dplyr::select(gene, ends_with("TP")) %>%
  filter(rowSums(dplyr::select(., ends_with("TP")), na.rm = TRUE) > 0)

cat("Number of genes expressed in both sexes in gonad:", length(mt_gonad$gene), "\n")
cat("Number of genes expressed in both sexes in tail:", length(mt_tail$gene), "\n")

# Pivot long
mt_gonad_long <- mt_gonad %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "potential_samp",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "gonad"
  )

mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "potential_samp",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

# Get effect size manually
get_log_rate_ratio <- function(df) {
  log(
    mean(df$potential_samp[df$sex == "male"]) /
      mean(df$potential_samp[df$sex == "female"])
  )
}

safe_ttest <- function(df) {
  tryCatch(
    t.test(potential_samp ~ sex, data = df)$p.value,
    error = function(e) NA_real_
  )
}

# Get sex effect size per gene
gonad_results_potential <- mt_gonad_long %>%
  group_by(gene) %>%
  summarise(
    log_rate_ratio = log2(
      mean(potential_samp[sex == "male"]) /
        mean(potential_samp[sex == "female"])
    ),
    p_value = safe_ttest(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    sig = !is.na(FDR) & FDR < 0.05
  )

glimpse(gonad_results_potential)

gonad_summary_potential <- gonad_results_potential %>%
  filter(is.finite(log_rate_ratio)) %>% # Remove anything sex-specific added back in with 0s
  summarise(
    mean_log = mean(log_rate_ratio, na.rm = TRUE),
    se_log = sd(log_rate_ratio, na.rm = TRUE) / sqrt(sum(!is.na(log_rate_ratio))),
    gene_count = sum(!is.na(log_rate_ratio))
  ) %>%
  mutate(
    mean_fold_change = exp(mean_log),
    lower_95 = exp(mean_log - 1.96 * se_log),
    upper_95 = exp(mean_log + 1.96 * se_log)
  )

volcano_df <- gonad_results_potential %>%
  filter(!is.na(p_value), is.finite(log_rate_ratio)) %>%
  mutate(
    category = case_when(
      FDR < 0.05 & log_rate_ratio > 0 ~ "More male isoforms",
      FDR < 0.05 & log_rate_ratio < 0 ~ "More female isoforms",
      TRUE ~ "Not significant"
    )
  )

# Mean location
mean_log_rate <- mean(volcano_df$log_rate_ratio, na.rm = TRUE)

# Get density height at mean position for placing the dot
dens <- density(volcano_df$log_rate_ratio, na.rm = TRUE)
mean_density_y <- approx(dens$x, dens$y, xout = mean_log_rate)$y

# Main volcano plot
p1 <- ggplot(volcano_df,
             aes(x = log_rate_ratio,
                 y = -log10(p_value),
                 color = category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "More male isoforms" = "#55BCC2",
      "More female isoforms" = "#E87D72",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = "log2FC (M:F)\nIsoform richness",
    y = expression(-log[10]("p-value"))
  ) +
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, 13)
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Density plot
p2 <- ggplot(volcano_df,
             aes(x = log_rate_ratio)) +
  geom_density(
    fill = "grey70",
    alpha = 0.5,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = mean_log_rate, y = 0),
    color = "red",
    size = 0.5
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

# Combine
p2 / p1 + plot_layout(heights = c(0.6, 1.7))

cat(
  "The mean effect size of sex (male/female) on potential in gonad is ",
  round(gonad_summary_potential$mean_fold_change, 3),
  " with 95% CI (",
  round(gonad_summary_potential$lower_95, 3), ", ",
  round(gonad_summary_potential$upper_95, 3), ").\n",
  sep = ""
)

### TAIL

# Get sex effect size per gene
tail_results_potential <- mt_tail_long %>%
  group_by(gene) %>%
  summarise(
    log_rate_ratio = log2(
      mean(potential_samp[sex == "male"]) /
        mean(potential_samp[sex == "female"])
    ),
    p_value = safe_ttest(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    sig = !is.na(FDR) & FDR < 0.05
  )

glimpse(tail_results_potential)

tail_summary_potential <- tail_results_potential %>%
  filter(is.finite(log_rate_ratio)) %>% # Remove anything sex-specific added back in with 0s
  summarise(
    mean_log = mean(log_rate_ratio, na.rm = TRUE),
    se_log = sd(log_rate_ratio, na.rm = TRUE) / sqrt(sum(!is.na(log_rate_ratio))),
    gene_count = sum(!is.na(log_rate_ratio))
  ) %>%
  mutate(
    mean_fold_change = exp(mean_log),
    lower_95 = exp(mean_log - 1.96 * se_log),
    upper_95 = exp(mean_log + 1.96 * se_log)
  )

volcano_df <- tail_results_potential %>%
  filter(!is.na(p_value), is.finite(log_rate_ratio)) %>%
  mutate(
    category = case_when(
      FDR < 0.05 & log_rate_ratio > 0 ~ "More male isoforms",
      FDR < 0.05 & log_rate_ratio < 0 ~ "More female isoforms",
      TRUE ~ "Not significant"
    )
  )

# Mean location
mean_log_rate <- mean(volcano_df$log_rate_ratio, na.rm = TRUE)

# Get density height at mean position for placing the dot
dens <- density(volcano_df$log_rate_ratio, na.rm = TRUE)
mean_density_y <- approx(dens$x, dens$y, xout = mean_log_rate)$y

# Main volcano plot
p1 <- ggplot(volcano_df,
             aes(x = log_rate_ratio,
                 y = -log10(p_value),
                 color = category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "More male isoforms" = "#55BCC2",
      "More female isoforms" = "#E87D72",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = "log2FC (M:F)\nIsoform richness",
    y = expression(-log[10]("p-value"))
  ) +
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, 6)
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Density plot
p2 <- ggplot(volcano_df,
             aes(x = log_rate_ratio)) +
  geom_density(
    fill = "grey70",
    alpha = 0.5,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = mean_log_rate, y = 0),
    color = "red",
    size = 0.5
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

# Combine
p2 / p1 + plot_layout(heights = c(0.6, 1.7))


cat(
  "The mean effect size of sex (male/female) on potential in tail is ",
  round(tail_summary_potential$mean_fold_change, 3),
  " with 95% CI (",
  round(tail_summary_potential$lower_95, 3), ", ",
  round(tail_summary_potential$upper_95, 3), ").\n",
  sep = ""
)

################################################################################
### (3) Leave one out analysis richness
################################################################################

glimpse(mt_samp)

cat("Number of genes expressed in both sexes in gonad:", length(mt_gonad$gene), "\n")
cat("Number of genes expressed in both sexes in tail:", length(mt_tail$gene), "\n")

# Pivot long
mt_gonad_long <- mt_gonad %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "potential_samp",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "gonad"
  )

mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "potential_samp",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

# Get effect size with a leave one out (loo) analysis
loo_effect <- function(df, tissue_name){
  
  samples <- unique(df$sample)
  
  purrr::map_dfr(samples, function(drop_sample){
    
    # Remove one replicate
    df_loo <- df %>%
      filter(sample != drop_sample)
    
    # Gene-level effect sizes
    results <- df_loo %>%
      group_by(gene) %>%
      summarise(
        log_rate_ratio = log2(
          mean(potential_samp[sex == "male"]) /
            mean(potential_samp[sex == "female"])
        ),
        .groups = "drop"
      ) %>%
      filter(is.finite(log_rate_ratio))
    
    # Overall mean effect
    summary <- results %>%
      summarise(
        mean_log = mean(log_rate_ratio),
        se_log = sd(log_rate_ratio)/sqrt(n()),
        gene_count = n()
      )
    
    tibble(
      tissue = tissue_name,
      omitted_sample = drop_sample,
      mean_log = summary$mean_log,
      se_log = summary$se_log,
      lower = summary$mean_log - 1.96 * summary$se_log,
      upper = summary$mean_log + 1.96 * summary$se_log,
      gene_count = summary$gene_count
    )
  })
}

gonad_loo <- loo_effect(mt_gonad_long, "gonad")
tail_loo <- loo_effect(mt_tail_long, "tail")

loo_results <- bind_rows(gonad_loo, tail_loo)

loo_results

full_results <- bind_rows(
  tibble(
    tissue = "gonad",
    omitted_sample = "None",
    mean_log = gonad_summary_potential$mean_log,
    se_log = gonad_summary_potential$se_log
  ),
  tibble(
    tissue = "tail",
    omitted_sample = "None",
    mean_log = tail_summary_potential$mean_log,
    se_log = tail_summary_potential$se_log
  )
) %>%
  mutate(
    lower = mean_log - 1.96 * se_log,
    upper = mean_log + 1.96 * se_log,
    type = "All data"
  )

plot_df <- bind_rows(
  loo_results %>%
    mutate(type = "Leave-one-out"),
  full_results %>%
    mutate(type = "All data")
)

# Rename tissues
plot_df <- plot_df %>%
  mutate(
    tissue = recode(
      tissue,
      "tail" = "Tail",
      "gonad" = "Gonad"
    )
  )

# Rename samples
sample_names <- c(
  "X1FTP" = "RFT1",
  "X2FTP" = "RFT2",
  "X3FTP" = "RFT3",
  "X4FTP" = "RFT4",
  "X5FTP" = "RFT5",
  "X6FTP" = "RFT6",
  "X1FGP" = "RFG1",
  "X2FGP" = "RFG2",
  "X3FGP" = "RFG3",
  "X4FGP" = "RFG4",
  "X5FGP" = "RFG5",
  "X6FGP" = "RFG6",
  "X1MTP" = "RMT1",
  "X2MTP" = "RMT2",
  "X3MTP" = "RMT3",
  "X4MTP" = "RMT4",
  "X5MTP" = "RMT5",
  "X6MTP" = "RMT6",
  "X1MGP" = "RMG1",
  "X2MGP" = "RMG2",
  "X3MGP" = "RMG3",
  "X4MGP" = "RMG4",
  "X5MGP" = "RMG5",
  "X6MGP" = "RMG6",
  "None"  = ""
)

plot_df <- plot_df %>%
  mutate(
    omitted_sample = recode(as.character(omitted_sample), !!!sample_names)
  )

# Put the full analysis ("") first, then samples in order
sample_order <- c(
  "",
  paste0("RFT", 1:6),
  paste0("RMT", 1:6),
  paste0("RFG", 1:6),
  paste0("RMG", 1:6)
)

plot_df <- plot_df %>%
  mutate(
    omitted_sample = factor(omitted_sample,
                            levels = rev(sample_order)),
    tissue = factor(tissue,
                    levels = c("Tail", "Gonad"))
  )

xlim <- max(abs(c(plot_df$lower, plot_df$upper)), na.rm = TRUE)

ggplot(plot_df,
       aes(
         y = omitted_sample,
         x = mean_log,
         colour = type
       )) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey60") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.15) +
  geom_point(size = 3) +
  facet_wrap(~tissue, scales = "free_y") +
  scale_colour_manual(values = c(
    "All data" = "red3",
    "Leave-one-out" = "black"
  )) +
  scale_x_continuous(limits = c(-xlim, xlim)) +
  labs(
    x = expression(log[2]("M/F Isoform richness")),
    y = "Replicate omitted",
    colour = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "right"
  )

################################################################################
### (4) Box plot for isoform richness (replicate means)
################################################################################

glimpse(mt_gonad_long)
glimpse(mt_tail_long)

mt_long <- bind_rows(mt_gonad_long, mt_tail_long)

# Already removed genes not expressed in both sexes and 0s back in

# Aggregate data to find the mean richness per replicate
replicate_means <- mt_long %>%
  # Flatten rows so each unique gene counts once per sample
  distinct(sample, tissue, sex, gene, .keep_all = TRUE) %>%
  filter(is.finite(potential_samp)) %>%
  group_by(sample, tissue, sex) %>% 
  summarise(
    mean_richness = mean(potential_samp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = case_when(
      tissue == "gonad" ~ "Gonad",
      tissue == "tail"  ~ "Tail",
      TRUE ~ tissue
    ),
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male")) # Clean factor levels for ordering
  )

# Create a combined factor group to compare all 4 conditions directly
test_data <- replicate_means %>%
  mutate(group = interaction(tissue, sex, sep = "_"))

# Run pairwise Wilcoxon tests between all group combinations
pairwise_wilcox <- pairwise.wilcox.test(
  x = test_data$mean_richness, 
  g = test_data$group, 
  p.adjust.method = "BH" # Benjamini-Hochberg adjustment for multiple comparisons
)

# Extract p-values and convert to a flat named vector format for multcompView
p_matrix <- pairwise_wilcox$p.value
p_vector <- na.omit(as.vector(p_matrix))

# Construct pairwise hyphenated comparison names (e.g., "Tail_female-Gonad_female")
names_vector <- c()
for(col in 1:ncol(p_matrix)) {
  for(row in 1:nrow(p_matrix)) {
    if(!is.na(p_matrix[row, col])) {
      names_vector <- c(names_vector, paste0(rownames(p_matrix)[row], "-", colnames(p_matrix)[col]))
    }
  }
}
names(p_vector) <- names_vector

# Generate compact letters display (shared letter = no sig diff)
cld_letters <- multcompLetters(p_vector, threshold = 0.05)$Letters

# Bind the letter flags back into a clean summary data frame for ggplot
letters_df <- data.frame(group = names(cld_letters), letter = cld_letters) %>%
  tidyr::separate(group, into = c("tissue", "sex"), sep = "_") %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

# Review the resulting letters layout in your console
print(letters_df)

# Manually define the letters for your plot to correct for order
letters_df <- data.frame(
  tissue = c("Tail", "Tail", "Gonad", "Gonad"),
  sex    = c("female", "male", "female", "male"),
  letter = c("a", "b", "c", "d") 
) %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

# Plot the sample-level means separated by sex, colored by group
richness_replicate_plot <- ggplot(replicate_means, aes(x = sex, y = mean_richness, color = sex)) +
  # Boxplot with fill, matching outline color, and hidden outliers
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.4, size = 0.7) +
  # Overlay individual replicate points colored dynamically by sex
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), 
             alpha = 0.8, size = 1, show.legend = FALSE) +
  # Overlay significance lettering based on Wilcoxon test matrix
  geom_text(data = letters_df, 
            aes(x = sex, y = 3.9, label = letter), # Positioned safely at the top view range
            size = 4, 
            color = "black",
            show.legend = FALSE) +
  # Facet by tissue cleanly side by side
  facet_wrap(~tissue, scales = "free_x") +
  scale_color_manual(
    values = c("female" = "#E87D72", "male" = "#55BCC2"),
    name = NULL
  ) +
  labs(
    y = "Isoform richness",
    x = NULL
  ) +
  scale_y_continuous(breaks = c(2, 3, 4)) +
  coord_cartesian(ylim = c(2.5, 4)) +
  theme_classic() + 
  theme(
    strip.background = element_rect(fill = 'grey90', color = NA),
    strip.text = element_text(),
    axis.text.x = element_blank(),  # Removes female/male text labels from the x-axis
    axis.ticks.x = element_blank(), # Removes x-axis ticks
    legend.position = "right",      
    legend.title = NULL)

richness_replicate_plot

################################################################################
### (5) Orf barplots and forest plots (isoform richness)
################################################################################

# We want 3 categories, and then subcategories within each category with similar hues:
# (1) Non-coding
# (2) Potential nonsense-mediated decay
#   (2a) Retained intron
#   (2b) Premature STOP
# (3) Coding
#   (3a) Unique ORF (unique AA sequences) = ORF variant
#   (X) Unique CDS (unique NT, redundant AA) = CDS variant # SKIP THIS ONE (only 10-20 isoforms)
#   (3b) Redundant CDS (redundant NT, only changed at the UTRs) = UTR variant

# Genes present in both sexes within a tissue
shared_genes <- mt_treat %>%
  distinct(tissue, sex, gene) %>%
  group_by(tissue, gene) %>%
  summarise(n_sexes = n_distinct(sex), .groups = "drop") %>%
  filter(n_sexes == 2)

# Assign types at the transcript level, keeping only shared genes
pot_codingpotential <- mt_treat %>%
  semi_join(shared_genes, by = c("tissue", "gene")) %>%
  dplyr::select(sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type))

# Instead of grouping by gene, we group by the plotting variables
summary_df <- pot_codingpotential %>%
  group_by(tissue, sex, iso_type) %>%
  summarise(total_count = n(), .groups = "drop_last") %>%
  mutate(
    prop = total_count / sum(total_count),
    label = paste0(round(prop * 100, 1), "%")
  ) %>%
  ungroup()

# Factor ordering for consistent visualization
summary_df <- summary_df %>%
  mutate(
    iso_type = factor(iso_type,
                      levels = c(
                        "Non-coding",
                        "Retained intron",
                        "Premature STOP",
                        "Unique ORF",
                        "UTR variant"
                      )),
    tissue = factor(tissue, levels = c("tail", "gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

### Stacked bar plots together with labels
# summary_df <- summary_df %>%
#   group_by(tissue, sex) %>%
#   mutate(
#     prop_within = total_count / sum(total_count),
#     label = paste0(
#       scales::comma(total_count),
#       " (",
#       round(prop_within * 100, 1),
#       "%)" # BOTH ABSOLUTE AND PROP LABELLED
#     )
#   ) %>%
#   ungroup()

summary_df <- summary_df %>%
  group_by(tissue, sex) %>%
  mutate(
    prop_within = total_count / sum(total_count),
    label = paste0(
      scales::comma(total_count)
    )
  ) %>%
  ungroup()

# Rearrange factor levels so "Unique ORF" is first (the bottom of the stack)
summary_df <- summary_df %>%
  mutate(iso_type = factor(iso_type, levels = c(
    "Non-coding",               # Top of the stack
    "Retained intron", 
    "Premature STOP", 
    "UTR variant", 
    "Unique ORF"        # Bottom of the stack
  )))

print(summary_df)

# Plot
pot_orf <- ggplot(summary_df,
                  aes(x = sex, y = prop, fill = iso_type)) +
  geom_col(position = "fill", width = 0.9) +
  #  geom_text(
  #    aes(label = label),
  #    position = position_fill(vjust = 0.5),
  #    size = 3
  #  ) +
  facet_wrap(~ tissue, nrow = 1,
             labeller = labeller(tissue = c("tail" = "Tail", "gonad" = "Gonad"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "Non-coding"              = "#D4C2C6",
      "Retained intron" = "#858893",
      "Premature STOP"   = "#BB6C87",
      "UTR variant"     = "#F8DCA7",
      "Unique ORF"   = "#97C7C0"
    )
  ) +
  labs(
    x = NULL,
    y = "Isoforms",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = 'grey90', color = NA),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
pot_orf

### DO REPLICATE LEVEL ANALYSIS

### mt_samp has only shared genes but no 0s added back
glimpse(mt_samp)

# Assign types at the transcript level, keeping only shared genes
pot_codingpotential <- mt_samp %>%
  dplyr::select(sample, sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sample, sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type))


# Calculate replicate-level isoform proportions
isoform_prop <- pot_codingpotential %>%
  count(sample, sex, tissue, iso_type) %>%
  group_by(sample, tissue) %>%
  mutate(prop = n / sum(n)) %>%
  dplyr::select(sample, sex, tissue, iso_type, prop) %>%
  pivot_wider(
    names_from = iso_type,
    values_from = prop
  )


# Fit replicate-aware models within each tissue and isoform class
results <- isoform_prop %>%
  pivot_longer(
    cols = -c(sample, sex, tissue),
    names_to = "iso_type",
    values_to = "prop"
  ) %>%
  group_by(tissue, iso_type) %>%
  do({
    model <- glm(prop ~ sex, data = .)
    tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sexmale") %>%
  group_by(tissue) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0 ~ "Higher in males",
      estimate < 0 ~ "Higher in females"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )

# Reorder factors
results <- results %>%
  mutate(
    tissue = factor(
      tissue,
      levels = c("tail", "gonad")
    )
  )

# Plot
ggplot(results,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ tissue) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Higher in females" = "#E87D72",
      "Higher in males" = "#55BCC2",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in isoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

################################################################################
### (6) Make log2FC for potential vs. bias (isoform richness) GONAD ONLY
################################################################################

# Combine data and set order (Tail then Gonad)
mt_combined_pot <- bind_rows(
  gonad_results_potential %>% mutate(tissue = "Gonad"),
  tail_results_potential %>% mutate(tissue = "Tail")
) %>% 
  filter(is.finite(log_rate_ratio)) %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    log2_fc = log_rate_ratio / log(2)  # convert ln to log2 fold change
  )


# Keep only Gonad data
mt_gonad_pot <- mt_combined_pot %>%
  filter(tissue == "Gonad") %>%
  left_join(
    combined_dge %>%
      mutate(
        tissue = case_when(
          tissue == "gonad" ~ "Gonad",
          tissue == "tail" ~ "Tail",
          TRUE ~ tissue
        )
      ) %>%
      dplyr::select(gene, tissue, dge.bias.direction, logFC),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      is.na(dge.bias.direction) ~ "unbiased",
      dge.bias.direction == "male-biased" ~ "male-biased",
      dge.bias.direction == "female-biased" ~ "female-biased",
      TRUE ~ "unbiased"
    ),
    sex_bias = factor(
      sex_bias,
      levels = c("female-biased", "unbiased", "male-biased")
    )
  )

# Summary table
means_table_gonad <- mt_gonad_pot %>%
  group_by(sex_bias) %>%
  summarise(
    mean_log2 = mean(log2_fc, na.rm = TRUE),
    se_log2 = sd(log2_fc, na.rm = TRUE) / sqrt(sum(!is.na(log2_fc))),
    n = sum(!is.na(log2_fc)),
    .groups = "drop"
  ) %>%
  mutate(
    lower_95_log2 = mean_log2 - 1.96 * se_log2,
    upper_95_log2 = mean_log2 + 1.96 * se_log2,
    mean_fc = 2^mean_log2,
    lower_95_fc = 2^lower_95_log2,
    upper_95_fc = 2^upper_95_log2,
    label = paste0(
      round(mean_log2, 3), 
      " [", round(lower_95_log2, 3), ", ",
      round(upper_95_log2, 3), "]\n",
      "n = ", n
    )
  )

# Kruskal-Wallis test
kw_gonad <- mt_gonad_pot %>%
  kruskal_test(log2_fc ~ sex_bias)

kw_gonad

# Dunn posthoc
dunn_gonad <- mt_gonad_pot %>%
  dunn_test(
    log2_fc ~ sex_bias,
    p.adjust.method = "BH"
  )
dunn_gonad

# Significance letters
groups <- sort(unique(c(dunn_gonad$group1, dunn_gonad$group2)))

pmat <- matrix(
  1,
  nrow = length(groups),
  ncol = length(groups),
  dimnames = list(groups, groups)
)

for(i in seq_len(nrow(dunn_gonad))) {
  g1 <- dunn_gonad$group1[i]
  g2 <- dunn_gonad$group2[i]
  p  <- dunn_gonad$p.adj[i]
  
  pmat[g1, g2] <- p
  pmat[g2, g1] <- p
}

diag(pmat) <- 0

cld_gonad <- multcompView::multcompLetters(pmat < 0.01)$Letters

cld_df <- tibble(
  sex_bias = names(cld_gonad),
  Letter = unname(cld_gonad)
)

# Position letters
label_df_gonad <- mt_gonad_pot %>%
  group_by(sex_bias) %>%
  summarise(
    y = 3.8,
    .groups = "drop"
  ) %>%
  left_join(cld_df, by = "sex_bias")

label_df_gonad <-  label_df_gonad %>%
  mutate(Letter = case_when(
    sex_bias == "female-biased" ~ "a",
    sex_bias == "unbiased" ~ "b",
    sex_bias == "male-biased" ~ "c"
  ))

# Plot
pot_bias_gonad <- ggplot(
  mt_gonad_pot,
  aes(x = sex_bias, y = log2_fc, fill = sex_bias)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_half_violin(
    side = "r",
    alpha = 0.7,
    color = NA,
    width = 1.3
  ) +
  geom_half_boxplot(
    side = "l",
    width = 0.6,          # Harmonized to match diversity plot
    outlier.shape = NA,
    fill = "white",
    aes(color = sex_bias)
  ) +
  geom_point(
    data = means_table_gonad,
    aes(x = sex_bias, y = mean_log2),
    color = "red",
    size = 1.5
  ) +
  geom_text(
    data = label_df_gonad,
    aes(x = sex_bias, y = y, label = Letter),
    size = 4,
    na.rm = TRUE
  ) +
  labs(
    y = "log2FC (M:F)\nIsoform richness",
    x = NULL,
    color = NULL
  ) +
  scale_fill_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    )
  ) +
  scale_color_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    )
  ) +
  guides(fill = "none") +
  coord_cartesian(
    ylim = c(-2.5, 3.8)   # Keeping specific custom y-limits
  ) +
  scale_x_discrete(       # Added to match spacing/width of diversity plot
    expand = expansion(mult = 0.5)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

pot_bias_gonad

################################################################################
### (7) Box plot for isoform richness vs bias GONAD ONLY
################################################################################

# mt_long is the df with sex-limited genes removed and 0s put back in
glimpse(mt_long)

# Filter mt_samp to Gonad expressed genes
mt_binned_samp_gonad <- mt_long %>%
  mutate(
    tissue = case_when(
      tissue == "gonad" ~ "Gonad",
      tissue == "tail"  ~ "Tail",
      TRUE ~ tissue
    )
  ) %>%
  filter(
    tissue == "Gonad"
  ) %>%
  distinct(sample, tissue, sex, gene, .keep_all = TRUE) %>%
  filter(is.finite(potential_samp))

# Join edgeR bias information
mt_binned_samp_gonad <- mt_binned_samp_gonad %>%
  left_join(
    combined_dge %>%
      mutate(
        tissue = case_when(
          tissue == "gonad" ~ "Gonad",
          tissue == "tail"  ~ "Tail",
          TRUE ~ tissue
        )
      ) %>%
      filter(tissue == "Gonad") %>%
      dplyr::select(
        gene,
        tissue,
        dge.bias.direction,
        logFC
      ),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      is.na(dge.bias.direction) ~ "unbiased",
      dge.bias.direction == "male-biased" ~ "male-biased",
      dge.bias.direction == "female-biased" ~ "female-biased",
      TRUE ~ "unbiased"
    ),
    sex_bias = factor(
      sex_bias,
      levels = c("female-biased", "unbiased", "male-biased")
    )
  )

# Define factor order
bias_levels <- c(
  "female-biased",
  "unbiased",
  "male-biased"
)

# Replicate aggregation: mean richness per sample per bias bin
replicate_binned_means_gonad <- mt_binned_samp_gonad %>%
  group_by(sample, tissue, sex, sex_bias) %>%
  summarise(
    mean_richness = mean(potential_samp),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = "Gonad"),
    sex_bias = factor(sex_bias, levels = bias_levels),
    sex = factor(sex, levels = c("female", "male"))
  )

# Statistical dataset (Gonad only)
statistical_data_gonad <- replicate_binned_means_gonad %>%
  mutate(
    clean_bias = stringr::str_replace_all(sex_bias, "-", "_"),
    group = clean_bias
  )

# Pairwise Wilcoxon tests between bias categories
pairwise_wilcox_gonad <- pairwise.wilcox.test(
  x = statistical_data_gonad$mean_richness,
  g = statistical_data_gonad$group,
  p.adjust.method = "BH"
)

pairwise_wilcox_gonad

# Convert lower triangle matrix to full symmetric matrix
p_matrix_gonad <- pairwise_wilcox_gonad$p.value

all_groups_gonad <- unique(statistical_data_gonad$group)

full_p_matrix_gonad <- matrix(
  1,
  nrow = length(all_groups_gonad),
  ncol = length(all_groups_gonad),
  dimnames = list(all_groups_gonad, all_groups_gonad)
)

for(col in colnames(p_matrix_gonad)) {
  for(row in rownames(p_matrix_gonad)) {
    if(!is.na(p_matrix_gonad[row, col])) {
      full_p_matrix_gonad[row, col] <- p_matrix_gonad[row, col]
      full_p_matrix_gonad[col, row] <- p_matrix_gonad[row, col]
    }
  }
}

# Compact letter display
cld_letters_gonad <- multcompLetters(
  full_p_matrix_gonad,
  threshold = 0.05
)$Letters

# Format letters for plotting
letters_df_gonad <- tibble::enframe(
  cld_letters_gonad,
  name = "group",
  value = "letter"
) %>%
  mutate(
    sex_bias = stringr::str_replace_all(group, "_", "-"),
    sex_bias = factor(sex_bias, levels = bias_levels)
  ) %>%
  dplyr::select(sex_bias, letter)


print(letters_df_gonad)

# Plot
richness_bias_plot_gonad <- ggplot(
  replicate_binned_means_gonad,
  aes(x = sex_bias, y = mean_richness)
) +
  geom_boxplot(
    aes(color = sex_bias),
    fill = "transparent",
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.4,
    size = 0.7
  ) +
  scale_color_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    ),
    name = NULL
  ) +
  ggnewscale::new_scale_color() +
  geom_point(
    aes(color = sex),
    position = position_jitter(width = 0.15),
    alpha = 0.8,
    size = 1.2,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "female" = "#E87D72",
      "male"   = "#55BCC2"
    ),
    name = NULL
  ) +
  geom_text(
    data = letters_df_gonad,
    aes(
      x = sex_bias,
      y = 5.8,
      label = letter
    ),
    color = "black",
    size = 4,
    show.legend = FALSE
  ) +
  labs(
    y = "Isoform richness",
    x = NULL
  ) +
  scale_y_continuous(
    breaks = c(1.5, 3.5, 5.5)
  ) +
  coord_cartesian(
    xlim = c(1, 3),
    ylim = c(1.5, 5.8)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

richness_bias_plot_gonad

################################################################################
### (8) Isoform composition female biased vs. male biased (gonad only)
################################################################################

### DO REPLICATE LEVEL ANALYSIS

### mt_samp has only shared genes but no 0s added back
glimpse(mt_samp)

# Assign types at the transcript level, keeping only shared genes
pot_codingpotential <- mt_samp %>%
  dplyr::select(sample, sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sample, sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type)) %>%  
  full_join(
    combined_dge %>%
      filter(tissue == "gonad") %>%
      dplyr::select(
        gene,
        tissue,
        dge.is.biased,
        logFC
      ),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      dge.is.biased == FALSE ~ "unbiased",
      logFC > 1              ~ "male-biased",
      logFC < -1             ~ "female-biased",
      TRUE                   ~ "unbiased"
    )
  )
  
# Calculate replicate-level isoform proportions within sex-bias class
isoform_prop <- pot_codingpotential %>%
  filter(tissue == "gonad",
         sex_bias %in% c("male-biased", "female-biased")) %>%
  count(sample, tissue, sex_bias, iso_type) %>%
  group_by(sample, tissue, sex_bias) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# Fit models
results <- isoform_prop %>%
  group_by(tissue, iso_type) %>%
  do({
    model <- lm(prop ~ sex_bias, data = .)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sex_biasmale-biased") %>%
  group_by(tissue) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0  ~ "Enriched in male-biased genes",
      TRUE          ~ "Enriched in female-biased genes"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )

# Reorder factors
results <- results %>%
  mutate(
    tissue = factor(tissue,
                    levels = "gonad",
                    labels = "Gonad")
  )

# Plot
ggplot(results,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ tissue) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Enriched in female-biased genes" = "#E05A4E",
      "Enriched in male-biased genes" = "#08519c",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in isoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

### FACET BY SEX

# Calculate replicate-level isoform proportions within sex-bias class and sex
isoform_prop_sex <- pot_codingpotential %>%
  filter(tissue == "gonad",
         sex_bias %in% c("male-biased", "female-biased")) %>%
  count(sample, sex, tissue, sex_bias, iso_type) %>%
  group_by(sample, sex, tissue, sex_bias) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# Fit models separately by sex
results_sex <- isoform_prop_sex %>%
  group_by(tissue, sex, iso_type) %>%
  do({
    model <- lm(prop ~ sex_bias, data = .)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sex_biasmale-biased") %>%
  group_by(tissue, sex) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0  ~ "Enriched in male-biased genes",
      TRUE          ~ "Enriched in female-biased genes"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )


# Reorder factors
results_sex <- results_sex %>%
  mutate(
    tissue = factor(tissue,
                    levels = "gonad",
                    labels = "Gonad")
  )


# Plot faceted by sex
ggplot(results_sex,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ sex) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Enriched in female-biased genes" = "#E05A4E",
      "Enriched in male-biased genes" = "#08519c",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in isoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

################################################################################
### (9) Calculate effect size of sex for perplexity (isoform diversity)
################################################################################

glimpse(mt_samp)

# Collapse to gene level and pivot wide, sub in 0 for NA
### We've already removed genes not expressed in both sexes
mt_ready <- mt_samp %>%
  dplyr::select(gene, sample, pGene_sample) %>%
  distinct() %>%
  pivot_wider(
    names_from = sample, 
    values_from = pGene_sample,
    values_fill = 0 
  )  

# Split tissues safely with na.rm = TRUE
mt_gonad <- mt_ready %>%
  dplyr::select(gene, ends_with("GP")) %>%
  filter(rowSums(dplyr::select(., ends_with("GP")), na.rm = TRUE) > 0)

mt_tail <- mt_ready %>%
  dplyr::select(gene, ends_with("TP")) %>%
  filter(rowSums(dplyr::select(., ends_with("TP")), na.rm = TRUE) > 0)

cat("Number of genes expressed in gonad:", length(mt_gonad$gene), "\n")
cat("Number of genes expressed in tail:", length(mt_tail$gene), "\n")

# Pivot long and remove NAs instantly
mt_gonad_long <- mt_gonad %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE # <-- Drops rows where pGene_sample is NA
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "gonad"
  )

mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE # <-- Drops rows where pGene_sample is NA
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

# Get effect size manually
get_log_rate_ratio <- function(df) {
  log(
    mean(df$pGene_sample[df$sex == "male"]) /
      mean(df$pGene_sample[df$sex == "female"])
  )
}

safe_ttest <- function(df) {
  tryCatch(
    t.test(pGene_sample ~ sex, data = df)$p.value,
    error = function(e) NA_real_
  )
}

# Get sex effect size per gene
gonad_results_perplexity <- mt_gonad_long %>%
  group_by(gene) %>%
  summarise(
    log_rate_ratio = log2(
      mean(pGene_sample[sex == "male"]) /
        mean(pGene_sample[sex == "female"])
    ),
    p_value = safe_ttest(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    sig = !is.na(FDR) & FDR < 0.05
  )

glimpse(gonad_results_perplexity)

gonad_summary_perplexity <- gonad_results_perplexity %>%
  filter(is.finite(log_rate_ratio)) %>%
  summarise(
    mean_log = mean(log_rate_ratio, na.rm = TRUE),
    se_log = sd(log_rate_ratio, na.rm = TRUE) / sqrt(sum(!is.na(log_rate_ratio))),
    gene_count = sum(!is.na(log_rate_ratio))
  ) %>%
  mutate(
    mean_fold_change = 2^mean_log,
    lower_95 = 2^(mean_log - 1.96 * se_log),
    upper_95 = 2^(mean_log + 1.96 * se_log)
  )

volcano_df <- gonad_results_perplexity %>%
  filter(!is.na(p_value), is.finite(log_rate_ratio)) %>%
  mutate(
    category = case_when(
      FDR < 0.05 & log_rate_ratio > 0 ~ "More male isoforms",
      FDR < 0.05 & log_rate_ratio < 0 ~ "More female isoforms",
      TRUE ~ "Not significant"
    )
  )

# Mean location
mean_log_rate <- mean(volcano_df$log_rate_ratio, na.rm = TRUE)

# Get density height at mean position for placing the dot
dens <- density(volcano_df$log_rate_ratio, na.rm = TRUE)
mean_density_y <- approx(dens$x, dens$y, xout = mean_log_rate)$y

# Main volcano plot
p1 <- ggplot(volcano_df,
             aes(x = log_rate_ratio,
                 y = -log10(p_value),
                 color = category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "More male isoforms" = "#55BCC2",
      "More female isoforms" = "#E87D72",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = "log2FC (M:F)\nIsoform diversity",
    y = expression(-log[10]("p-value"))
  ) +
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, 13)
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Density plot
p2 <- ggplot(volcano_df,
             aes(x = log_rate_ratio)) +
  geom_density(
    fill = "grey70",
    alpha = 0.5,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = mean_log_rate, y = 0),
    color = "red",
    size = 0.5
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

# Combine
p2 / p1 + plot_layout(heights = c(0.6, 1.7))


cat(
  "The mean effect size of sex (male/female) on perplexity in gonad is ",
  round(gonad_summary_perplexity$mean_fold_change, 3),
  " with 95% CI (",
  round(gonad_summary_perplexity$lower_95, 3), ", ",
  round(gonad_summary_perplexity$upper_95, 3), ").\n",
  sep = ""
)

### Tail

# Pivot long and remove NAs instantly
mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE # <-- Drops rows where pGene_sample is NA
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE # <-- Drops rows where pGene_sample is NA
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

# Get effect size manually
get_log_rate_ratio <- function(df) {
  log(
    mean(df$pGene_sample[df$sex == "male"]) /
      mean(df$pGene_sample[df$sex == "female"])
  )
}

safe_ttest <- function(df) {
  tryCatch(
    t.test(pGene_sample ~ sex, data = df)$p.value,
    error = function(e) NA_real_
  )
}

# Get sex effect size per gene
tail_results_perplexity <- mt_tail_long %>%
  group_by(gene) %>%
  summarise(
    log_rate_ratio = log2(
      mean(pGene_sample[sex == "male"]) /
        mean(pGene_sample[sex == "female"])
    ),
    p_value = safe_ttest(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    sig = !is.na(FDR) & FDR < 0.05
  )

glimpse(tail_results_perplexity)

tail_summary_perplexity <- tail_results_perplexity %>%
  filter(is.finite(log_rate_ratio)) %>%
  summarise(
    mean_log = mean(log_rate_ratio, na.rm = TRUE),
    se_log = sd(log_rate_ratio, na.rm = TRUE) / sqrt(sum(!is.na(log_rate_ratio))),
    gene_count = sum(!is.na(log_rate_ratio))
  ) %>%
  mutate(
    mean_fold_change = 2^mean_log,
    lower_95 = 2^(mean_log - 1.96 * se_log),
    upper_95 = 2^(mean_log + 1.96 * se_log)
  )

volcano_df <- tail_results_perplexity %>%
  filter(!is.na(p_value), is.finite(log_rate_ratio)) %>%
  mutate(
    category = case_when(
      FDR < 0.05 & log_rate_ratio > 0 ~ "More male isoforms",
      FDR < 0.05 & log_rate_ratio < 0 ~ "More female isoforms",
      TRUE ~ "Not significant"
    )
  )

# Mean location
mean_log_rate <- mean(volcano_df$log_rate_ratio, na.rm = TRUE)

# Get density height at mean position for placing the dot
dens <- density(volcano_df$log_rate_ratio, na.rm = TRUE)
mean_density_y <- approx(dens$x, dens$y, xout = mean_log_rate)$y

# Main volcano plot
p1 <- ggplot(volcano_df,
             aes(x = log_rate_ratio,
                 y = -log10(p_value),
                 color = category)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c(
      "More male isoforms" = "#55BCC2",
      "More female isoforms" = "#E87D72",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = "log2FC (M:F)\nIsoform diversity",
    y = expression(-log[10]("p-value"))
  ) +
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, 13)
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Density plot
p2 <- ggplot(volcano_df,
             aes(x = log_rate_ratio)) +
  geom_density(
    fill = "grey70",
    alpha = 0.5,
    adjust = 1.5
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = mean_log_rate, y = 0),
    color = "red",
    size = 0.5
  ) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

# Combine
p2 / p1 + plot_layout(heights = c(0.6, 1.7))


cat(
  "The mean effect size of sex (male/female) on perplexity in tail is ",
  round(tail_summary_perplexity$mean_fold_change, 3),
  " with 95% CI (",
  round(tail_summary_perplexity$lower_95, 3), ", ",
  round(tail_summary_perplexity$upper_95, 3), ").\n",
  sep = ""
)

################################################################################
### (10) Leave one out analysis diversity
################################################################################

# Load sample level data
mt_samp <- read.csv(here::here("data",
                               "isoseq",
                               "11_mastertable",
                               "mastertable_sample.csv"))

glimpse(mt_samp)

### FILTER NOTICE: Remove any genes not expressed in both sexes
mt_samp <- mt_samp %>%
  group_by(tissue, gene) %>%
  filter(
    all(c("female", "male") %in%
          sex[geneTPM > 0])
  ) %>%
  ungroup()

### PUT BACK IN 0s FOR REPLICATES WITHOUT EXP
# Collapse to gene level and pivot wide, sub in 0 for NA
mt_ready <- mt_samp %>%
  dplyr::select(gene, sample, pGene_sample) %>%
  distinct() %>%
  pivot_wider(
    names_from = sample, 
    values_from = pGene_sample,
    values_fill = 0 
  )

# Split tissues safely with na.rm = TRUE
mt_gonad <- mt_ready %>%
  dplyr::select(gene, ends_with("GP")) %>%
  filter(rowSums(dplyr::select(., ends_with("GP")), na.rm = TRUE) > 0)

mt_tail <- mt_ready %>%
  dplyr::select(gene, ends_with("TP")) %>%
  filter(rowSums(dplyr::select(., ends_with("TP")), na.rm = TRUE) > 0)

cat("Number of genes expressed in both sexes in gonad:", length(mt_gonad$gene), "\n")
cat("Number of genes expressed in both sexes in tail:", length(mt_tail$gene), "\n")

# Pivot long
mt_gonad_long <- mt_gonad %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "gonad"
  )

mt_tail_long <- mt_tail %>%
  pivot_longer(
    cols = -gene, 
    names_to = "sample", 
    values_to = "pGene_sample",
    values_drop_na = TRUE 
  ) %>%
  mutate(
    sex = ifelse(str_detect(sample, "F"), "female", "male"),
    tissue = "tail"
  )

# Get effect size with a leave one out (loo) analysis
loo_effect <- function(df, tissue_name){
  
  samples <- unique(df$sample)
  
  purrr::map_dfr(samples, function(drop_sample){
    
    # Remove one replicate
    df_loo <- df %>%
      filter(sample != drop_sample)
    
    # Gene-level effect sizes
    results <- df_loo %>%
      group_by(gene) %>%
      summarise(
        log_rate_ratio = log2(
          mean(pGene_sample[sex == "male"]) /
            mean(pGene_sample[sex == "female"])
        ),
        .groups = "drop"
      ) %>%
      filter(is.finite(log_rate_ratio))
    
    # Overall mean effect
    summary <- results %>%
      summarise(
        mean_log = mean(log_rate_ratio),
        se_log = sd(log_rate_ratio)/sqrt(n()),
        gene_count = n()
      )
    
    tibble(
      tissue = tissue_name,
      omitted_sample = drop_sample,
      mean_log = summary$mean_log,
      se_log = summary$se_log,
      lower = summary$mean_log - 1.96 * summary$se_log,
      upper = summary$mean_log + 1.96 * summary$se_log,
      gene_count = summary$gene_count
    )
  })
}

gonad_loo <- loo_effect(mt_gonad_long, "gonad")
tail_loo <- loo_effect(mt_tail_long, "tail")

loo_results <- bind_rows(gonad_loo, tail_loo)

loo_results

full_results <- bind_rows(
  tibble(
    tissue = "gonad",
    omitted_sample = "None",
    mean_log = gonad_summary_perplexity$mean_log,
    se_log = gonad_summary_perplexity$se_log
  ),
  tibble(
    tissue = "tail",
    omitted_sample = "None",
    mean_log = tail_summary_perplexity$mean_log,
    se_log = tail_summary_perplexity$se_log
  )
) %>%
  mutate(
    lower = mean_log - 1.96 * se_log,
    upper = mean_log + 1.96 * se_log,
    type = "All data"
  )

plot_df <- bind_rows(
  loo_results %>%
    mutate(type = "Leave-one-out"),
  full_results %>%
    mutate(type = "All data")
)

# Rename tissues
plot_df <- plot_df %>%
  mutate(
    tissue = recode(
      tissue,
      "tail" = "Tail",
      "gonad" = "Gonad"
    )
  )

# Rename samples
sample_names <- c(
  "X1FTP" = "RFT1",
  "X2FTP" = "RFT2",
  "X3FTP" = "RFT3",
  "X4FTP" = "RFT4",
  "X5FTP" = "RFT5",
  "X6FTP" = "RFT6",
  "X1FGP" = "RFG1",
  "X2FGP" = "RFG2",
  "X3FGP" = "RFG3",
  "X4FGP" = "RFG4",
  "X5FGP" = "RFG5",
  "X6FGP" = "RFG6",
  "X1MTP" = "RMT1",
  "X2MTP" = "RMT2",
  "X3MTP" = "RMT3",
  "X4MTP" = "RMT4",
  "X5MTP" = "RMT5",
  "X6MTP" = "RMT6",
  "X1MGP" = "RMG1",
  "X2MGP" = "RMG2",
  "X3MGP" = "RMG3",
  "X4MGP" = "RMG4",
  "X5MGP" = "RMG5",
  "X6MGP" = "RMG6",
  "None"  = ""
)

plot_df <- plot_df %>%
  mutate(
    omitted_sample = recode(as.character(omitted_sample), !!!sample_names)
  )

# Put the full analysis ("") first, then samples in order
sample_order <- c(
  "",
  paste0("RFT", 1:6),
  paste0("RMT", 1:6),
  paste0("RFG", 1:6),
  paste0("RMG", 1:6)
)

plot_df <- plot_df %>%
  mutate(
    omitted_sample = factor(omitted_sample,
                            levels = rev(sample_order)),
    tissue = factor(tissue,
                    levels = c("Tail", "Gonad"))
  )

xlim <- max(abs(c(plot_df$lower, plot_df$upper)), na.rm = TRUE)

ggplot(plot_df,
       aes(
         y = omitted_sample,
         x = mean_log,
         colour = type
       )) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey60") +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.15) +
  geom_point(size = 3) +
  facet_wrap(~tissue, scales = "free_y") +
  scale_colour_manual(values = c(
    "All data" = "red3",
    "Leave-one-out" = "black"
  )) +
  scale_x_continuous(limits = c(-xlim, xlim)) +
  labs(
    x = expression(log[2]("M/F Isoform diversity")),
    y = "Replicate omitted",
    colour = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "right"
  )

################################################################################
### (11) Box plot for isoform diversity (replicate means)
################################################################################

### This is the df with sex-limited genes removed and 0s added back
mt_long <- bind_rows(mt_gonad_long, mt_tail_long)

# Aggregate data to find the mean diversity per replicate
replicate_means <- mt_long %>%
  semi_join(shared_genes, by = c("tissue", "gene")) %>% 
  # Flatten rows so each unique gene counts once per sample
  distinct(sample, tissue, sex, gene, .keep_all = TRUE) %>%
  filter(is.finite(pGene_sample)) %>%
  group_by(sample, tissue, sex) %>% 
  summarise(
    mean_diversity = mean(pGene_sample, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = case_when(
      tissue == "gonad" ~ "Gonad",
      tissue == "tail"  ~ "Tail",
      TRUE ~ tissue
    ),
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male")) # Clean factor levels for ordering
  )

# Create a combined factor group to compare all 4 conditions directly
test_data <- replicate_means %>%
  mutate(group = interaction(tissue, sex, sep = "_"))

# Run pairwise Wilcoxon tests between all group combinations
pairwise_wilcox <- pairwise.wilcox.test(
  x = test_data$mean_diversity, 
  g = test_data$group, 
  p.adjust.method = "BH" # Benjamini-Hochberg adjustment for multiple comparisons
)

# Extract p-values and convert to a flat named vector format for multcompView
p_matrix <- pairwise_wilcox$p.value
p_vector <- na.omit(as.vector(p_matrix))

# Construct pairwise hyphenated comparison names (e.g., "Tail_female-Gonad_female")
names_vector <- c()
for(col in 1:ncol(p_matrix)) {
  for(row in 1:nrow(p_matrix)) {
    if(!is.na(p_matrix[row, col])) {
      names_vector <- c(names_vector, paste0(rownames(p_matrix)[row], "-", colnames(p_matrix)[col]))
    }
  }
}
names(p_vector) <- names_vector

# Generate compact letters display (shared letter = no sig diff)
cld_letters <- multcompLetters(p_vector, threshold = 0.05)$Letters

# Bind the letter flags back into a clean summary data frame for ggplot
letters_df <- data.frame(group = names(cld_letters), letter = cld_letters) %>%
  tidyr::separate(group, into = c("tissue", "sex"), sep = "_") %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

# Review the resulting letters layout in your console
print(letters_df)

# Manually define the letters for your plot to correct for order
letters_df <- data.frame(
  tissue = c("Tail", "Tail", "Gonad", "Gonad"),
  sex    = c("female", "male", "female", "male"),
  letter = c("a", "b", "a", "c") 
) %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

# Plot the sample-level means separated by sex, colored by group
diversity_replicate_plot <- ggplot(replicate_means, aes(x = sex, y = mean_diversity, color = sex)) +
  # Boxplot with fill, matching outline color, and hidden outliers
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.4, size = 0.7) +
  # Overlay individual replicate points colored dynamically by sex
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), 
             alpha = 0.8, size = 1, show.legend = FALSE) +
  # Overlay significance lettering based on Wilcoxon test matrix
  geom_text(data = letters_df, 
            aes(x = sex, y = 1.97, label = letter), # Positioned safely at the top view range
            size = 4, 
            color = "black",
            show.legend = FALSE) +
  # Facet by tissue cleanly side by side
  facet_wrap(~tissue, scales = "free_x") +
  scale_color_manual(
    values = c("female" = "#E87D72", "male" = "#55BCC2"),
    name = NULL
  ) +
  labs(
    y = "Isoform diversity",
    x = NULL
  ) +
  scale_y_continuous(breaks = c(1.6, 2)) +
  coord_cartesian(ylim = c(1.6, 2)) +
  theme_classic() + 
  theme(
    strip.background = element_rect(fill = 'grey90', color = NA),
    strip.text = element_text(),
    axis.text.x = element_blank(),  # Removes female/male text labels from the x-axis
    axis.ticks.x = element_blank(), # Removes x-axis ticks
    legend.position = "right",      
    legend.title = NULL)

diversity_replicate_plot

################################################################################
### (12) Orf barplots (isoform diversity)
################################################################################

# We want 3 categories, and then subcategories within each category with similar hues:
# (1) Non-coding
# (2) Potential nonsense-mediated decay
#   (2a) Retained intron
#   (2b) Premature STOP
# (3) Coding
#   (3a) Unique ORF (unique AA sequences) = ORF variant
#   (X) Unique CDS (unique NT, redundant AA) = CDS variant # SKIP THIS ONE (only 10-20 isoforms)
#   (3b) Redundant CDS (redundant NT, only changed at the UTRs) = UTR variant

# Genes present in both sexes within a tissue
shared_genes <- mt_treat %>%
  distinct(tissue, sex, gene) %>%
  group_by(tissue, gene) %>%
  summarise(n_sexes = n_distinct(sex), .groups = "drop") %>%
  filter(n_sexes == 2)

# Assign types at the transcript level, keeping only shared genes
perp_codingpotential <- mt_treat %>%
  semi_join(shared_genes, by = c("tissue", "gene")) %>%
  filter(effGene_treat == T) %>% ### Effective isoforms only
  dplyr::select(sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type))

# Instead of grouping by gene, we group by the plotting variables
summary_df <- perp_codingpotential %>%
  group_by(tissue, sex, iso_type) %>%
  summarise(total_count = n(), .groups = "drop_last") %>%
  mutate(
    prop = total_count / sum(total_count),
    label = paste0(round(prop * 100, 1), "%")
  ) %>%
  ungroup()

# Factor ordering for consistent visualization
summary_df <- summary_df %>%
  mutate(
    iso_type = factor(iso_type,
                      levels = c(
                        "Non-coding",
                        "Retained intron",
                        "Premature STOP",
                        "Unique ORF",
                        "UTR variant"
                      )),
    tissue = factor(tissue, levels = c("tail", "gonad")),
    sex    = factor(sex, levels = c("female", "male"))
  )

### Stacked bar plots together with labels
# summary_df <- summary_df %>%
#   group_by(tissue, sex) %>%
#   mutate(
#     prop_within = total_count / sum(total_count),
#     label = paste0(
#       scales::comma(total_count),
#       " (",
#       round(prop_within * 100, 1),
#       "%)" # BOTH ABSOLUTE AND PROP LABELLED
#     )
#   ) %>%
#   ungroup()

summary_df <- summary_df %>%
  group_by(tissue, sex) %>%
  mutate(
    prop_within = total_count / sum(total_count),
    label = paste0(
      scales::comma(total_count)
    )
  ) %>%
  ungroup()

# Rearrange factor levels so "Unique ORF" is first (the bottom of the stack)
summary_df <- summary_df %>%
  mutate(iso_type = factor(iso_type, levels = c(
    "Non-coding",               # Top of the stack
    "Retained intron", 
    "Premature STOP", 
    "UTR variant", 
    "Unique ORF"        # Bottom of the stack
  )))

print(summary_df)

# Plot
perp_orf <- ggplot(summary_df,
                  aes(x = sex, y = prop, fill = iso_type)) +
  geom_col(position = "fill", width = 0.9) +
  #  geom_text(
  #    aes(label = label),
  #    position = position_fill(vjust = 0.5),
  #    size = 3
  #  ) +
  facet_wrap(~ tissue, nrow = 1,
             labeller = labeller(tissue = c("tail" = "Tail", "gonad" = "Gonad"))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "Non-coding"              = "#D4C2C6",
      "Retained intron" = "#858893",
      "Premature STOP"   = "#BB6C87",
      "UTR variant"     = "#F8DCA7",
      "Unique ORF"   = "#97C7C0"
    )
  ) +
  labs(
    x = NULL,
    y = "Isoforms",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = 'grey90', color = NA),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
perp_orf

### DO REPLICATE LEVEL ANALYSIS

### mt_samp has only shared genes but no 0s added back
glimpse(mt_samp)

# Assign transcript-level isoform types using shared genes
perp_codingpotential <- mt_samp %>%
  filter(effGene_sample == TRUE) %>%
  dplyr::select(sample, sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sample, sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type))


# Calculate replicate-level isoform proportions
isoform_prop <- perp_codingpotential %>%
  count(sample, sex, tissue, iso_type) %>%
  group_by(sample, tissue) %>%
  mutate(prop = n / sum(n)) %>%
  dplyr::select(sample, sex, tissue, iso_type, prop) %>%
  pivot_wider(
    names_from = iso_type,
    values_from = prop
  )


# Fit replicate-aware models within each tissue and isoform class
results <- isoform_prop %>%
  pivot_longer(
    cols = -c(sample, sex, tissue),
    names_to = "iso_type",
    values_to = "prop"
  ) %>%
  group_by(tissue, iso_type) %>%
  do({
    model <- glm(prop ~ sex, data = .)
    tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sexmale") %>%
  group_by(tissue) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0 ~ "Higher in males",
      estimate < 0 ~ "Higher in females"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )

# Reorder factors
results <- results %>%
  mutate(
    tissue = factor(
      tissue,
      levels = c("tail", "gonad")
    )
  )

# Plot
ggplot(results,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ tissue) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Higher in females" = "#E87D72",
      "Higher in males" = "#55BCC2",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in effective isoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

################################################################################
### (13) Make log2FC for perplexity vs. bias (isoform diversity)
################################################################################

# Combine data and set order (Tail then Gonad)
mt_combined_perp <- bind_rows(
  gonad_results_perplexity %>% mutate(tissue = "Gonad"),
  tail_results_perplexity %>% mutate(tissue = "Tail")
) %>% 
  filter(is.finite(log_rate_ratio)) %>%
  mutate(
    tissue = factor(tissue, levels = c("Tail", "Gonad")),
    log2_fc = log_rate_ratio / log(2)  # convert ln to log2 fold change
  )

# Keep only Gonad data
mt_gonad_perp <- mt_combined_perp %>%
  filter(tissue == "Gonad") %>%
  left_join(
    combined_dge %>%
      mutate(tissue = case_when(
        tissue == "gonad" ~ "Gonad",
        tissue == "tail" ~ "Tail",
        TRUE ~ tissue
      )) %>%
      dplyr::select(gene, tissue, dge.is.biased, logFC),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      dge.is.biased == FALSE ~ "unbiased",
      logFC > 1              ~ "male-biased",
      logFC < -1             ~ "female-biased",
      TRUE                   ~ "unbiased"
    ),
    sex_bias = factor(
      sex_bias,
      levels = c("female-biased", "unbiased", "male-biased")
    )
  )

# Summary table
means_table_gonad_perp <- mt_gonad_perp %>%
  group_by(sex_bias) %>%
  summarise(
    mean_log2 = mean(log2_fc, na.rm = TRUE),
    se_log2 = sd(log2_fc, na.rm = TRUE) / sqrt(sum(!is.na(log2_fc))),
    n = sum(!is.na(log2_fc)),
    .groups = "drop"
  ) %>%
  mutate(
    lower_95_log2 = mean_log2 - 1.96 * se_log2,
    upper_95_log2 = mean_log2 + 1.96 * se_log2,
    mean_fc = 2^mean_log2,
    lower_95_fc = 2^lower_95_log2,
    upper_95_fc = 2^upper_95_log2,
    label = paste0(
      round(mean_log2, 3),
      " [", round(lower_95_log2, 3), ", ",
      round(upper_95_log2, 3), "]\n",
      "n = ", n
    )
  )

# Kruskal-Wallis test
kw_gonad_perp <- mt_gonad_perp %>%
  kruskal_test(log2_fc ~ sex_bias)

kw_gonad_perp

# Dunn posthoc test
dunn_gonad_perp <- mt_gonad_perp %>%
  dunn_test(
    log2_fc ~ sex_bias,
    p.adjust.method = "BH"
  )

dunn_gonad_perp

# Generate CLD letters
groups <- sort(unique(c(
  dunn_gonad_perp$group1,
  dunn_gonad_perp$group2
)))

pmat <- matrix(
  1,
  nrow = length(groups),
  ncol = length(groups),
  dimnames = list(groups, groups)
)

for(i in seq_len(nrow(dunn_gonad_perp))) {
  g1 <- dunn_gonad_perp$group1[i]
  g2 <- dunn_gonad_perp$group2[i]
  p  <- dunn_gonad_perp$p.adj[i]
  
  pmat[g1, g2] <- p
  pmat[g2, g1] <- p
}

diag(pmat) <- 0

cld_gonad_perp <- multcompView::multcompLetters(pmat < 0.01)$Letters

cld_df_gonad_perp <- tibble::tibble(
  sex_bias = names(cld_gonad_perp),
  Letter = unname(cld_gonad_perp)
)

# Position significance letters
label_df_gonad_perp <- mt_gonad_perp %>%
  group_by(sex_bias) %>%
  summarise(
    y = 2.5,   # move letters near the top of the plot
    .groups = "drop"
  ) %>%
  left_join(cld_df_gonad_perp, by = "sex_bias") %>%
  mutate(Letter = na_if(Letter, ""))

print(label_df_gonad_perp)

label_df_gonad_perp <- label_df_gonad_perp %>%
  mutate(Letter = case_when(
    sex_bias == "female-biased" ~ "a",
    sex_bias == "unbiased" ~ "b",
    sex_bias == "male-biased" ~ "c"
  ))

# Plot
perp_bias_gonad <- ggplot(
  mt_gonad_perp,
  aes(x = sex_bias, y = log2_fc, fill = sex_bias)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_half_violin(
    side = "r",
    alpha = 0.7,
    color = NA,
    width = 1.3
  ) +
  geom_half_boxplot(
    side = "l",
    width = 0.6,          # Harmonized to match richness plot
    outlier.shape = NA,
    fill = "white",
    aes(color = sex_bias)
  ) +
  geom_point(
    data = means_table_gonad_perp,
    aes(x = sex_bias, y = mean_log2),
    color = "red",
    size = 1.5
  ) +
  geom_text(
    data = label_df_gonad_perp,
    aes(x = sex_bias, y = y, label = Letter),
    size = 4,
    na.rm = TRUE
  ) +
  labs(
    y = "log2FC (M:F)\nIsoform diversity",
    x = NULL,
    color = NULL
  ) +
  scale_fill_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    )
  ) +
  scale_color_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    )
  ) +
  guides(fill = "none") +
  coord_cartesian(
    ylim = c(-1.2, 2.5)   # Keeping specific custom y-limits
  ) +
  scale_x_discrete(       # Keeping matching spacing/width
    expand = expansion(mult = 0.5)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

perp_bias_gonad

pot_bias_gonad <- pot_bias_gonad +
  theme(legend.position = "none")

perp_bias_gonad <- perp_bias_gonad +
  theme(legend.position = "none")

### COMBINED

combined_bias_gonad <- pot_bias_gonad + perp_bias_gonad +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

combined_bias_gonad

################################################################################
### (14) Box plot for isoform diversity vs bias GONAD ONLY
################################################################################

# mt_long is the df with sex-limited genes removed and 0s put back in
glimpse(mt_long)

# Filter mt_samp to Gonad expressed genes
mt_binned_samp_gonad <- mt_long %>%
  mutate(
    tissue = case_when(
      tissue == "gonad" ~ "Gonad",
      tissue == "tail"  ~ "Tail",
      TRUE ~ tissue
    )
  ) %>%
  filter(
    tissue == "Gonad"
  ) %>%
  distinct(sample, tissue, sex, gene, .keep_all = TRUE) %>%
  filter(is.finite(pGene_sample))

# Join Gonad edgeR bias calls
mt_binned_samp_gonad <- mt_binned_samp_gonad %>%
  left_join(
    combined_dge %>%
      mutate(
        tissue = case_when(
          tissue == "gonad" ~ "Gonad",
          tissue == "tail"  ~ "Tail",
          TRUE ~ tissue
        )
      ) %>%
      filter(tissue == "Gonad") %>%
      dplyr::select(
        gene,
        tissue,
        dge.is.biased,
        logFC
      ),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      dge.is.biased == FALSE ~ "unbiased",
      logFC > 1              ~ "male-biased",
      logFC < -1             ~ "female-biased",
      TRUE                   ~ "unbiased"
    )
  )

# Aggregate by gene
mt_binned_samp_gonad <- mt_binned_samp_gonad %>%
  dplyr::select(gene, sex, tissue, sex, sex_bias, sample, pGene_sample)

# Fill in 0s
sample_info <- mt_binned_samp_gonad %>%
  distinct(sample, sex, tissue)

replicate_complete <- mt_binned_samp_gonad %>%
  dplyr::select(sample, gene, sex_bias, pGene_sample) %>%   # <-- no sex or tissue
  distinct() %>%
  pivot_wider(
    names_from = sample,
    values_from = pGene_sample
  ) %>%
  mutate(across(-c(gene, sex_bias), ~replace_na(.x, 0))) %>%
  pivot_longer(
    -c(gene, sex_bias),
    names_to = "sample",
    values_to = "pGene_sample"
  ) %>%
  left_join(sample_info, by = "sample")

# Factor levels
bias_levels <- c(
  "female-biased",
  "unbiased",
  "male-biased"
)

# Replicate aggregation
replicate_binned_means_gonad <- replicate_complete %>%
  group_by(sample, tissue, sex, sex_bias) %>%
  summarise(
    mean_diversity = mean(pGene_sample, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = "Gonad"),
    sex_bias = factor(sex_bias, levels = bias_levels),
    sex = factor(sex, levels = c("female", "male"))
  )

# Statistical testing
statistical_data_gonad <- replicate_binned_means_gonad %>%
  mutate(
    clean_bias = stringr::str_replace_all(sex_bias, "-", "_"),
    group = clean_bias
  )

pairwise_wilcox_gonad <- pairwise.wilcox.test(
  x = statistical_data_gonad$mean_diversity,
  g = statistical_data_gonad$group,
  p.adjust.method = "BH"
)

pairwise_wilcox_gonad

# Convert p-values to full symmetric matrix
p_matrix_gonad <- pairwise_wilcox_gonad$p.value

all_groups_gonad <- unique(statistical_data_gonad$group)

full_p_matrix_gonad <- matrix(
  1,
  nrow = length(all_groups_gonad),
  ncol = length(all_groups_gonad),
  dimnames = list(all_groups_gonad, all_groups_gonad)
)

for(col in colnames(p_matrix_gonad)) {
  for(row in rownames(p_matrix_gonad)) {
    if(!is.na(p_matrix_gonad[row, col])) {
      full_p_matrix_gonad[row, col] <- p_matrix_gonad[row, col]
      full_p_matrix_gonad[col, row] <- p_matrix_gonad[row, col]
    }
  }
}

# Compact letters
cld_letters_gonad <- multcompLetters(
  full_p_matrix_gonad,
  threshold = 0.05
)$Letters

letters_df_gonad <- tibble::enframe(
  cld_letters_gonad,
  name = "group",
  value = "letter"
) %>%
  mutate(
    sex_bias = stringr::str_replace_all(group, "_", "-"),
    sex_bias = factor(sex_bias, levels = bias_levels)
  ) %>%
  dplyr::select(
    sex_bias,
    letter
  )

print(letters_df_gonad)

letters_df_gonad <- letters_df_gonad %>%
  mutate(
    letter = case_when(
      sex_bias == "female-biased" ~ "a",
      sex_bias == "unbiased" ~ "b",
      sex_bias == "male-biased" ~ "c"
    )
  )

# Plot with points colored by sex
diversity_bias_plot_gonad <- ggplot(
  replicate_binned_means_gonad,
  aes(
    x = sex_bias,
    y = mean_diversity
  )
) +
  # Boxplot outlines colored by bias category
  geom_boxplot(
    aes(color = sex_bias),
    fill = "transparent",
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.4,
    size = 0.7
  ) +
  scale_color_manual(
    values = c(
      "male-biased"   = "#08519c",
      "unbiased"      = "#bdbdbd",
      "female-biased" = "#E05A4E"
    ),
    name = NULL
  ) +
  # New color scale for sex points
  ggnewscale::new_scale_color() +
  geom_point(
    aes(color = sex),
    position = position_jitter(width = 0.15),
    alpha = 0.8,
    size = 1.2,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "female" = "#E87D72",
      "male"   = "#55BCC2"
    ),
    name = NULL
  ) +
  geom_text(
    data = letters_df_gonad,
    aes(
      x = sex_bias,
      y = 2.3,
      label = letter
    ),
    color = "black",
    size = 4,
    show.legend = FALSE
  ) +
  labs(
    y = "Isoform diversity",
    x = NULL
  ) +
  coord_cartesian(
    ylim = c(
      1.5,
      2.3
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

diversity_bias_plot_gonad


combined_bias_gonad <- richness_bias_plot_gonad + diversity_bias_plot_gonad+
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

combined_bias_gonad

################################################################################
### (15) Eff isoform composition female biased vs. male biased (gonad only)
################################################################################

### DO REPLICATE LEVEL ANALYSIS

### mt_samp has only shared genes but no 0s added back
glimpse(mt_samp)

# Assign types at the transcript level, keeping only shared genes
perp_codingpotential <- mt_samp %>%
  filter(effGene_sample == TRUE) %>%
  dplyr::select(sample, sex, tissue, gene, transcript, ptc_status, orf_identity) %>%
  group_by(sample, sex, tissue, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"         ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type)) %>%  
  full_join(
    combined_dge %>%
      filter(tissue == "gonad") %>%
      dplyr::select(
        gene,
        tissue,
        dge.is.biased,
        logFC
      ),
    by = c("gene", "tissue")
  ) %>%
  mutate(
    sex_bias = case_when(
      dge.is.biased == FALSE ~ "unbiased",
      logFC > 1              ~ "male-biased",
      logFC < -1             ~ "female-biased",
      TRUE                   ~ "unbiased"
    )
  )

# Calculate replicate-level isoform proportions within sex-bias class
isoform_prop <- perp_codingpotential %>%
  filter(tissue == "gonad",
         sex_bias %in% c("male-biased", "female-biased")) %>%
  count(sample, tissue, sex_bias, iso_type) %>%
  group_by(sample, tissue, sex_bias) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# Fit models
results <- isoform_prop %>%
  group_by(tissue, iso_type) %>%
  do({
    model <- lm(prop ~ sex_bias, data = .)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sex_biasmale-biased") %>%
  group_by(tissue) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0  ~ "Enriched in male-biased genes",
      TRUE          ~ "Enriched in female-biased genes"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )

# Reorder factors
results <- results %>%
  mutate(
    tissue = factor(tissue,
                    levels = "gonad",
                    labels = "Gonad")
  )

# Plot
ggplot(results,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ tissue) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Enriched in female-biased genes" = "#E05A4E",
      "Enriched in male-biased genes" = "#08519c",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in effective \nisoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

### Do it again with sex faceted

# Calculate replicate-level isoform proportions within sex-bias class and sex
isoform_prop_sex <- perp_codingpotential %>%
  filter(tissue == "gonad",
         sex_bias %in% c("male-biased", "female-biased")) %>%
  count(sample, sex, tissue, sex_bias, iso_type) %>%
  group_by(sample, sex, tissue, sex_bias) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


# Fit models separately by sex
results_sex <- isoform_prop_sex %>%
  group_by(tissue, sex, iso_type) %>%
  do({
    model <- lm(prop ~ sex_bias, data = .)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term == "sex_biasmale-biased") %>%
  group_by(tissue, sex) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      p_adj >= 0.05 ~ "Not significant",
      estimate > 0  ~ "Enriched in male-biased genes",
      TRUE          ~ "Enriched in female-biased genes"
    ),
    iso_type = factor(
      iso_type,
      levels = c(
        "Non-coding",
        "Retained intron",
        "Premature STOP",
        "UTR variant",
        "Unique ORF"
      )
    )
  )


# Reorder factors
results_sex <- results_sex %>%
  mutate(
    tissue = factor(tissue,
                    levels = "gonad",
                    labels = "Gonad"),
    sex = factor(sex)
  )


# Plot faceted by sex
ggplot(results_sex,
       aes(x = iso_type,
           y = estimate,
           color = direction)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ),
    width = 0.2
  ) +
  facet_wrap(~ sex) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c(
      "Enriched in female-biased genes" = "#E05A4E",
      "Enriched in male-biased genes" = "#08519c",
      "Not significant" = "grey70"
    )
  ) +
  labs(
    x = NULL,
    y = "Change in effective\nisoform composition",
    color = NULL
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text = element_text(colour = "black"),
    
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.ticks.x = element_blank(),
    
    legend.position = "none"
  )

################################################################################
### (16) Bias orfs (richness and diversity)
################################################################################

### Coding potential - isoforms

# First, extract the gene-level bias classifications specifically for Gonad
gonad_bias_map <- mt_treat %>%
  filter(tissue %in% c("gonad", "Gonad")) %>%
  distinct(gene, dge.is.biased, logFC) %>%
  mutate(
    sex_bias = case_when(
      dge.is.biased == FALSE ~ "unbiased",
      logFC > 1              ~ "male-biased",
      logFC < -1             ~ "female-biased",
      TRUE                   ~ "unbiased"
    ),
    sex_bias = factor(sex_bias, levels = bias_levels)
  )

# Build the dataset, ensuring we process orf_first within sex before collapsing
summary_df_bias <- mt_treat %>%
  semi_join(shared_genes, by = c("tissue", "gene")) %>%
  filter(tissue %in% c("gonad", "Gonad")) %>%
  inner_join(gonad_bias_map, by = "gene") %>%
  dplyr::select(sex, gene, transcript, ptc_status, orf_identity, sex_bias) %>%
  
  # CRITICAL: Keep sex here so orf_first is correctly calculated!
  group_by(sex, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"          ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type)) %>%
  
  # NOW collapse the sexes together by grouping purely on bias class and iso_type
  group_by(sex_bias, iso_type) %>%
  summarise(total_count = n(), .groups = "drop_last") %>%
  mutate(
    prop = total_count / sum(total_count)
  ) %>%
  ungroup() %>%
  mutate(
    iso_type = factor(iso_type, levels = c(
      "Non-coding",       
      "Retained intron", 
      "Premature STOP", 
      "UTR variant", 
      "Unique ORF"        
    ))
  )

pot_orf_bias <- ggplot(summary_df_bias, aes(x = sex_bias, y = prop, fill = iso_type)) +
  geom_col(position = "fill", width = 0.6) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Non-coding"      = "#D4C2C6",
      "Retained intron" = "#858893",
      "Premature STOP"  = "#BB6C87",
      "UTR variant"     = "#F8DCA7",
      "Unique ORF"      = "#97C7C0"
    )
  ) +
  labs(
    x = NULL, 
    y = "Isoforms", 
    fill = "Isoform Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

print(pot_orf_bias)

### Coding potential - effective isoforms

# Build the dataset, adding the effective isoform filter stage
summary_df_bias <- mt_treat %>%
  semi_join(shared_genes, by = c("tissue", "gene")) %>%
  filter(tissue %in% c("gonad", "Gonad")) %>%
  inner_join(gonad_bias_map, by = "gene") %>%
  
  # NEW: Filter specifically for effective isoforms
  filter(effGene_treat == TRUE) %>% 
  
  dplyr::select(sex, gene, transcript, ptc_status, orf_identity, sex_bias) %>%
  
  # Group by sex and gene to identify uniquely translated transcripts correctly
  group_by(sex, gene) %>%
  arrange(orf_identity, transcript, .by_group = TRUE) %>%
  mutate(orf_first = !duplicated(orf_identity)) %>%
  ungroup() %>%
  
  mutate(
    iso_type = case_when(
      ptc_status == "non-coding"          ~ "Non-coding",
      ptc_status == "intron_retention"   ~ "Retained intron",
      ptc_status == "nmd"                ~ "Premature STOP",
      ptc_status == "coding" & orf_first ~ "Unique ORF",
      ptc_status == "coding" & !orf_first ~ "UTR variant",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(iso_type)) %>%
  # Collapse the data by bias class and structural form
  group_by(sex_bias, iso_type) %>%
  summarise(total_count = n(), .groups = "drop_last") %>%
  mutate(
    prop = total_count / sum(total_count)
  ) %>%
  ungroup() %>%
  mutate(
    iso_type = factor(iso_type, levels = c(
      "Non-coding",       
      "Retained intron", 
      "Premature STOP", 
      "UTR variant", 
      "Unique ORF"        
    ))
  )

perp_orf_bias <- ggplot(summary_df_bias, aes(x = sex_bias, y = prop, fill = iso_type)) +
  geom_col(position = "fill", width = 0.6) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), 
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Non-coding"      = "#D4C2C6",
      "Retained intron" = "#858893",
      "Premature STOP"  = "#BB6C87",
      "UTR variant"     = "#F8DCA7",
      "Unique ORF"      = "#97C7C0"
    )
  ) +
  labs(
    x = NULL, 
    y = "Effective Isoforms", 
    fill = "Isoform Type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

print(perp_orf_bias)




