# Run edgeR for DGE analysis

# Only running on genes expressed in both sexes

# (1) Load files
# (2) Create design matrix
# (3) Differential Gene Expression (DGE)
# (4) Differential Transcript Expression (DTE)
# (5) Compare DGE and DTE
# (6) Differential Transcript Usage (aka. differential splicing; DTU)
# (7) Venn diagram at gene level for each tissue
# (8) Venn diagram at transcript level for each tissue

library(edgeR)
library(dplyr)
library(ggVennDiagram)
library(VennDiagram)
library(tidyverse)
library(ggpointdensity)
library(scales)
library(dunn.test)
library(ggsignif)
library(brms)
library(cmdstanr)
library(limma)
library(ggvenn)
library(PNWColors)
library(vcd)

################################################################################
### (1) Load files
################################################################################

# Load edgeR objects
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "01_eR_genegonad.RData"))
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "02_eR_genetail.RData"))
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "01_eR_trangonad.RData"))
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "02_eR_trantail.RData"))

# Load master table
source("load_mt.R")

# Get isoforms from genes expressed in both sexes
tran_both_sexes <- mt_samp %>%
  filter(geneTPM > 0) %>%
  group_by(tissue, transcript, gene, sex) %>%
  summarise(
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  group_by(tissue, transcript, gene) %>%
  summarise(
    n_sexes = n_distinct(sex),
    .groups = "drop"
  ) %>%
  filter(n_sexes == 2)

# Split transcript IDs by tissue
trans_gonad <- tran_both_sexes %>%
  filter(tissue == "gonad") %>%
  pull(transcript) %>%
  unique()
genes_gonad <- tran_both_sexes %>%
  filter(tissue == "gonad") %>%
  pull(gene) %>%
  unique()
trans_tail <- tran_both_sexes %>%
  filter(tissue == "tail") %>%
  pull(transcript) %>%
  unique()
genes_tail <- tran_both_sexes %>%
  filter(tissue == "tail") %>%
  pull(gene) %>%
  unique()

# Filter dge objects
dge_gene_gonad <- dge_gene_gonad[
  rownames(dge_gene_gonad) %in% genes_gonad,
]
dge_gene_tail <- dge_gene_tail[
  rownames(dge_gene_tail) %in% genes_tail,
]
dge_tran_gonad <- dge_tran_gonad[
  rownames(dge_tran_gonad) %in% trans_gonad,
]
dge_tran_tail <- dge_tran_tail[
  rownames(dge_tran_tail) %in% trans_tail,
]

# Also load metadata for batch effects
# Contains: cell, sample_id, sex, tissue, selection, read depth (from flnc.bam)
metadata <- read.csv(here::here("data",
                                "isoseq",
                                "metadata.csv"))

# Your vector of all sample names
all_samples <- c("1FGP", "1FTP", "1MGP", "1MTP",
                 "2FGP", "2FTP", "2MGP", "2MTP",
                 "3FGP", "3FTP", "3MGP", "3MTP",
                 "4FGP", "4FTP", "4MGP", "4MTP",
                 "5FGP", "5FTP", "5MGP", "5MTP",
                 "6FGP", "6FTP", "6MGP", "6MTP")

# Combine into a metadata data frame
sample_info <- data.frame(
  sample = all_samples,
  sex = sex <- factor(ifelse(grepl("F", all_samples), "Female", "Male")),
  tissue = tissue <- factor(ifelse(grepl("T", all_samples), "Tail", "Gonad")),
  treatment = factor(paste0(sex, "_", tissue))
)

# Add batch effects from metadata (cell) and selection column
sample_info <- sample_info %>%
  left_join(metadata %>% 
              dplyr::select(sample_id, selection, cell),
            by = c("sample" = "sample_id"))

# Clean up 
rm(all_samples)
rm(metadata)
rm(sex)
rm(tissue)
rm(mt_all)
rm(mt_tissue)
rm(mt_treat)

################################################################################
### (2) Create design matrix
################################################################################

# Subset sample_info into gonad and gonad
sample_info_gonad <- subset(sample_info, tissue == "Gonad")
sample_info_gonad <- droplevels(sample_info_gonad)

sample_info_tail <- subset(sample_info, tissue == "Tail")
sample_info_tail <- droplevels(sample_info_tail)

# Design matrix
# Add batch effects in order to test them within model
design_g <- model.matrix(~ 0 + sex, data = sample_info_gonad)
colnames(design_g) <- c(
  "Female", 
  "Male"
)
rownames(design_g) <- sample_info_gonad$sample
design_t <- model.matrix(~ 0 + sex, data = sample_info_tail)
colnames(design_t) <- c(
  "Female", 
  "Male"
)

# Make contrast matrix
contrast_matrix_t <- makeContrasts(
  Male_vs_Female = Male - Female,
  levels = colnames(design_t)
)
contrast_matrix_g <- makeContrasts(
  Male_vs_Female = Male - Female,
  levels = colnames(design_g)
)

# Clean up
rm(sample_info_gonad)
rm(sample_info_tail)
rm(sample_info)

################################################################################
### (3) Differential Gene Expression (DGE)
################################################################################

# Sanity check
dim(dge_gene_gonad)
dim(dge_gene_tail)

# Estimate disp
dge_gene_gonad <- estimateDisp(dge_gene_gonad, design_g)
dge_gene_tail <- estimateDisp(dge_gene_tail, design_t)

# Fit model - here we use glmQLFit 
fit_g <- glmQLFit(dge_gene_gonad, design_g)
fit_t <- glmQLFit(dge_gene_tail, design_t)

# Run tests
gonad_test_gene <- glmTreat(fit_g, 
                            contrast=contrast_matrix_g[, "Male_vs_Female"],
                            lfc = 1)
tail_test_gene <- glmTreat(fit_t, 
                           contrast=contrast_matrix_t[, "Male_vs_Female"],
                           lfc = 1)

# Take a look
summary(decideTests(gonad_test_gene))
# 1*Female -1*Male
# Down               3567
# NotSig             5938
# Up                 3703
summary(decideTests(tail_test_gene))
# 1*Female -1*Male
# Down                  8
# NotSig            12480
# Up                  195

# Visualize
plotMD(gonad_test_gene)
abline(h = c(-1, 1), col = "blue")
plotMD(tail_test_gene)
abline(h = c(-1, 1), col = "blue")

# Extract results

# p-value is the raw sig score, FDR is after BH correction
# so FDR<0.05 is significant 
dge_gonad_results <- edgeR::topTags(gonad_test_gene, n = Inf)$table
dge_tail_results <- edgeR::topTags(tail_test_gene, n = Inf)$table

# Annotate bias
# Define significance based on both FDR and effect size threshold
dge_gonad_results$dge.is.biased <- ifelse(
  dge_gonad_results$FDR < 0.05 & abs(dge_gonad_results$logFC) >= 1, T, F)
dge_tail_results$dge.is.biased <- ifelse(
  dge_tail_results$FDR < 0.05 & abs(dge_tail_results$logFC) >= 1, T, F)

# Bias direction only for significant genes (with effect size threshold)
dge_gonad_results$dge.bias.direction <- "unbiased"  # initialize
dge_gonad_results$dge.bias.direction[dge_gonad_results$dge.is.biased == T & dge_gonad_results$logFC > 1] <- "male-biased"
dge_gonad_results$dge.bias.direction[dge_gonad_results$dge.is.biased == T & dge_gonad_results$logFC < -1] <- "female-biased"
dge_tail_results$dge.bias.direction <- "unbiased"  
dge_tail_results$dge.bias.direction[dge_tail_results$dge.is.biased == T & dge_tail_results$logFC > 1] <- "male-biased"
dge_tail_results$dge.bias.direction[dge_tail_results$dge.is.biased == T & dge_tail_results$logFC < -1] <- "female-biased"

# I don't think Overdispersion means anything here, it's just leftover from
# the transcript collapse into gene
dge_gonad_results$EffectiveLength <- NULL
dge_gonad_results$Overdispersion <- NULL
dge_tail_results$EffectiveLength <- NULL
dge_tail_results$Overdispersion <- NULL

# Clean up
rm(dge_gene_gonad)
rm(dge_gene_tail)
rm(gonad_test_gene)
rm(tail_test_gene)

# Save filtered DGE results
write.csv(dge_gonad_results,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "01_dge_gonad_long.csv"),
          row.names = F)
write.csv(dge_tail_results,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "02_dge_tail_long.csv"),
          row.names = F)

################################################################################
### (4) Differential Transcript Expression (DTE)
################################################################################

# Estimate dispersions with the quasi-likelihood pipeline

# Sanity check
dim(dge_tran_gonad)
dim(dge_tran_tail)

# Plot BCV (biological coefficient of variation - squareroot of dispersion) 
# - if you were to plot for the unscaled
# counts (before accounting for overdispersion) you will see the y-axis is higher
bcv_g <- estimateDisp(dge_tran_gonad, design_g)
plotBCV(bcv_g)
bcv_t <- estimateDisp(dge_tran_tail, design_t)
plotBCV(bcv_t)

# "Note that this NB dispersion estimation step (estimateDisp) is now optional as all the NB dispersion estimates
# will not be used further under the latest quasi-likelihood (QL) pipeline.
# For DTE analyses, we still recommend using the quasi-likelihood (QL) pipeline for stricter
# error rate control by accounting for the uncertainty associated with the dispersion estimation.
# To this end, we estimate QL dispersions via glmQLFit and visualize the estimated values with
# plotQLDisp." - edgeR manual

# Estimate dispersions with glmQLFit 
# (Genewise Negative Binomial Generalized Linear Models with Quasi-Dispersion Estimation)
# glmQLFit is similar to glmFit except that it also estimates a quasi-dispersion 
# for each transcript. It calls the limma function squeezeVar to conduct empirical 
# Bayes moderation of the transcript-wise quasi-dispersions. If robust=TRUE, 
# then the robust hyperparameter estimation features of squeezeVar are used (Phipson et al, 2016).

# "Under the latest QL pipeline, technical and biological variation can be estimated simultaneously 
# in the glmQLFit function. This returns a DGEGLM object containing the estimated values
# of the GLM coefficients for each gene, an overall estimate of the squared biological coefficient
# of variation, as well as the fitted mean-QL dispersion trend, the squeezed QL estimates and
# the prior degrees of freedom (df). These can be visualized with the plotQLDisp function"

# This fits the negative binomial GLM to your data, using the design matrix 
# (which includes the batch variable). You must run this step before any hypothesis testing.
fit_g <- glmQLFit(dge_tran_gonad, 
                  design_g,
                  robust = T) # setting robust = T is strongly recommended in edgeR doc
fit_t <- glmQLFit(dge_tran_tail, 
                  design_t,
                  robust = T) 

# Plot the mean-QL dispersion trend and the individual gene QL dispersion estimates
plotQLDisp(fit_g)
plotQLDisp(fit_t)

# Test for differential transcript expression (DTE)

# Use glmQLFTest to conduct genewise quasi F-tests for a given coefficient or coefficient contrast,
# glmQLFTest() tests: H₀: logFC = 0
# glmTreat() tests:   H₀: logFC ≤ lfc, e.g. logFC ≤ 1 - this is equivalent to 2x expression
# So it's more stringent — useful when you care about biologically meaningful fold changes
gonad_test_transcript <- glmTreat(fit_g, 
                                  contrast = contrast_matrix_g[,"Male_vs_Female"],
                                  lfc = 1)
tail_test_transcript <- glmTreat(fit_t, 
                                 contrast = contrast_matrix_t[,"Male_vs_Female"],
                                 lfc = 1)

# The top set of most significant transcripts can be examined with topTags. Here, a positive log-fold
# change represents genes that are up in FemaleGonad/gonad over MaleGonad/gonad. Multiplicity correction is
# performed by applying the Benjamini-Hochberg method on the p-values, to control the false
# discovery rate (FDR)
# topTags(gonad_test_transcript)
# topTags(tail_test_transcript)

# The total number of DE genes in each direction at a FDR of 5% can be examined with decideTests.
# Uses BH method to correct.
summary(decideTests(gonad_test_transcript))
summary(decideTests(tail_test_transcript))

# The differential expression test results can be visualized using an MD plot. The log-fold
# change for each gene is plotted against the average abundance, i.e., logTPM in the result table
# above. Significantly DE genes at a FDR of 5% are highlighted.
plotMD(gonad_test_transcript)
abline(h=c(-1,1), col="blue")
plotMD(tail_test_transcript)
abline(h=c(-1,1), col="blue")

# Extract full DTE results 
dte_gonad_results <- gonad_test_transcript$table
dte_tail_results <- tail_test_transcript$table

# Add BH correction
dte_gonad_results$FDR <-p.adjust(dte_gonad_results$PValue, method = "BH")
dte_tail_results$FDR <- p.adjust(dte_tail_results$PValue, method = "BH")

# Extract transcript and gene info
transcript_info_g <- gonad_test_transcript$genes[, c("TranscriptID", "GeneID")]
transcript_info_t <- tail_test_transcript$genes[, c("TranscriptID", "GeneID")]

# Combine into one data frame by row order (they should be aligned)
dte_gonad_results <- cbind(transcript_info_g, dte_gonad_results)
dte_tail_results <- cbind(transcript_info_t, dte_tail_results)

# Add is biased column (FDR < 0.05 and absolute logFC >= 1)
dte_gonad_results$dte.is.biased <- ifelse(dte_gonad_results$FDR < 0.05 
                                          & abs(dte_gonad_results$logFC) >= 1, T, F)
dte_tail_results$dte.is.biased <- ifelse(dte_tail_results$FDR < 0.05 
                                         & abs(dte_tail_results$logFC) >= 1, T, F)


# Bias direction only for significant genes (with effect size threshold)
dte_gonad_results$dte.bias.direction <- "unbiased"
dte_gonad_results$dte.bias.direction[dte_gonad_results$dte.is.biased == T & dte_gonad_results$logFC > 1] <- "male-biased"
dte_gonad_results$dte.bias.direction[dte_gonad_results$dte.is.biased == T & dte_gonad_results$logFC < -1] <- "female-biased"

dte_tail_results$dte.bias.direction <- "unbiased"  # initialize with NA
dte_tail_results$dte.bias.direction[dte_tail_results$dte.is.biased == T & dte_tail_results$logFC > 1] <- "male-biased"
dte_tail_results$dte.bias.direction[dte_tail_results$dte.is.biased == T & dte_tail_results$logFC < -1] <- "female-biased"

# Clean up
rm(gonad_test_transcript)
rm(tail_test_transcript)

# Save DTE results
write.csv(dte_gonad_results,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "01_dte_gonad_long.csv"),
          row.names = F)
write.csv(dte_tail_results,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "02_dte_tail_long.csv"),
          row.names = F)

################################################################################
### (5) Compare DGE and DTE
################################################################################

### First build contingency tables

# TAIL

# Join Transcript (DIE) and Gene (DGE) data for Tail
tail_paired <- dte_tail_results %>%
  select(TranscriptID, GeneID, dte.bias.direction) %>%
  inner_join(
    dge_tail_results %>% select(GeneID, dge.bias.direction),
    by = "GeneID"
  ) %>%
  mutate(
    dge.bias.direction = factor(dge.bias.direction, levels = c("female-biased", "unbiased", "male-biased")),
    dte.bias.direction = factor(dte.bias.direction, levels = c("female-biased", "unbiased", "male-biased"))
  )

# Build contingency table (DGE = rows, DIE = columns)
tail_matrix <- table(
  DGE = tail_paired$dge.bias.direction, 
  DIE = tail_paired$dte.bias.direction
)

# Calculate Row Percentages (margin = 1 for Gene Perspective)
tail_row_pct <- prop.table(tail_matrix, margin = 1) * 100
tail_totals  <- rowSums(tail_matrix)

# Format cells as: Count (Row %)
tail_formatted_cells <- matrix(
  paste0(tail_matrix, " (", sprintf("%.2f", tail_row_pct), "%)"),
  nrow = nrow(tail_matrix),
  dimnames = dimnames(tail_matrix)
)

# Combine with total counts column
tail_final_df <- as.data.frame(tail_formatted_cells) %>%
  mutate(`Total isoform count` = tail_totals)

cat("--- TAIL: GENE PERSPECTIVE TABLE (ROW PERCENTAGES) ---\n")
print(tail_final_df)

cat("\n--- COHEN'S KAPPA TEST (TAIL) ---\n")
tail_kappa <- vcd::Kappa(tail_matrix)
print(tail_kappa)

# GONAD

# Join Transcript (DIE) and Gene (DGE) data for Gonad
gonad_paired <- dte_gonad_results %>%
  select(TranscriptID, GeneID, dte.bias.direction) %>%
  inner_join(
    dge_gonad_results %>% select(GeneID, dge.bias.direction),
    by = "GeneID"
  ) %>%
  mutate(
    dge.bias.direction = factor(dge.bias.direction, levels = c("female-biased", "unbiased", "male-biased")),
    dte.bias.direction = factor(dte.bias.direction, levels = c("female-biased", "unbiased", "male-biased"))
  )

# Build contingency table (DGE = rows, DIE = columns)
gonad_matrix <- table(
  DGE = gonad_paired$dge.bias.direction, 
  DIE = gonad_paired$dte.bias.direction
)

# Calculate Row Percentages (margin = 1 for Gene Perspective)
gonad_row_pct <- prop.table(gonad_matrix, margin = 1) * 100
gonad_totals  <- rowSums(gonad_matrix)

# Format cells as: Count (Row %)
gonad_formatted_cells <- matrix(
  paste0(gonad_matrix, " (", sprintf("%.2f", gonad_row_pct), "%)"),
  nrow = nrow(gonad_matrix),
  dimnames = dimnames(gonad_matrix)
)

# Combine with total counts column
gonad_final_df <- as.data.frame(gonad_formatted_cells) %>%
  mutate(`Total isoform count` = gonad_totals)

cat("\n--- GONAD: GENE PERSPECTIVE TABLE (ROW PERCENTAGES) ---\n")
print(gonad_final_df)

cat("\n--- COHEN'S KAPPA TEST (GONAD) ---\n")
gonad_kappa <- vcd::Kappa(gonad_matrix)
print(gonad_kappa)

### Prepare for plotting

# Prepare DTE data (Transcript level)
plot_dte <- bind_rows(
  dte_gonad_results %>% mutate(Tissue = "Gonad"),
  dte_tail_results %>% mutate(Tissue = "Tail")
) %>%
  dplyr::group_by(Tissue, dte.bias.direction) %>%
  tally() %>%
  dplyr::rename(Bias = dte.bias.direction) %>%
  mutate(Level = "Transcript")

# Prepare DGE data (Gene level)
plot_dge <- bind_rows(
  dge_gonad_results %>% mutate(Tissue = "Gonad"),
  dge_tail_results %>% mutate(Tissue = "Tail")
) %>%
  group_by(Tissue, dge.bias.direction) %>%
  tally() %>%
  dplyr::rename(Bias = dge.bias.direction) %>%
  mutate(Level = "Gene")

# Combine and set Factor Levels
final_plot_data <- bind_rows(plot_dte, plot_dge) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Tail", "Gonad")),
    Bias = factor(Bias, levels = c("female-biased", "unbiased", "male-biased")),
    Level = factor(Level, levels = c("Gene", "Transcript"))
  )

# Relabel levels
final_plot_data <- final_plot_data %>%
  dplyr::mutate(
    Level = dplyr::recode(Level,
                          "Gene" = "Gene",
                          "Transcript" = "Isoform")
  )

# Gene plot
gene_df <- final_plot_data %>% filter(Level == "Gene")

p_gene <- ggplot(gene_df, aes(x = Tissue, y = n, fill = Bias)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.85),
    width = 0.85,
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = scales::comma(n)),
    position = position_dodge(width = 0.85),
    vjust = -0.5,
    size = 3.2,
    color = "#2c3e50"
  ) +
  facet_wrap(~Tissue, scales = "free_x", labeller = as_labeller(str_to_title)) +
  scale_fill_manual(values = c(
    "female-biased" = "#E05A4E",
    "male-biased"   = "#08519c",
    "unbiased"      = "#bdbdbd"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = NULL,
    y = "Genes", 
    title = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "plain", size = 11, color = "black"),
    # FIX: Completely removed the bottom labels and ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.5)
  )

# Isoform plot
iso_df <- final_plot_data %>% filter(Level == "Isoform")

p_isoform <- ggplot(iso_df, aes(x = Tissue, y = n, fill = Bias)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.85),
    width = 0.85,
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = scales::comma(n)),
    position = position_dodge(width = 0.85),
    vjust = -0.5,
    size = 3.2,
    color = "#2c3e50"
  ) +
  facet_wrap(~Tissue, scales = "free_x", labeller = as_labeller(str_to_title)) +
  scale_fill_manual(values = c(
    "female-biased" = "#E05A4E",
    "male-biased"   = "#08519c",
    "unbiased"      = "#bdbdbd"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = NULL,
    y = "Isoforms", 
    title = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "plain", size = 11, color = "black"),
    # FIX: Completely removed the bottom labels and ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey95", linewidth = 0.5),
    legend.position = "none"
  )

p_gene
p_isoform

### Try presenting it together

# Prepare Gene-level classifications
gene_class_gonad <- dge_gonad_results %>% 
  dplyr::select(GeneID, Gene_Bias = dge.bias.direction)

gene_class_tail <- dge_tail_results %>% 
  dplyr::select(GeneID, Gene_Bias = dge.bias.direction)

# Join Gene classifications to the Transcript data
# This assigns a "Gene_Bias" label to every individual transcript
dte_nested_gonad <- dte_gonad_results %>%
  inner_join(gene_class_gonad, by = "GeneID") %>%
  mutate(Tissue = "Gonad")

dte_nested_tail <- dte_tail_results %>%
  inner_join(gene_class_tail, by = "GeneID") %>%
  mutate(Tissue = "Tail")

# Combine and Plot
plot_nested_data <- bind_rows(dte_nested_gonad, dte_nested_tail) %>%
  mutate(
    Tissue = factor(Tissue, levels = c("Tail", "Gonad")),
    Gene_Bias = factor(Gene_Bias, levels = c("female-biased", "unbiased", "male-biased")),
    dte.bias.direction = factor(dte.bias.direction, levels = c("female-biased", "unbiased", "male-biased"))
  )

# Calculate proportions
nested_summary <- plot_nested_data %>%
  count(Tissue, Gene_Bias, dte.bias.direction) %>%
  group_by(Tissue, Gene_Bias) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

# Plot proportion plot
nested_props <- ggplot(
  nested_summary,
  aes(
    x = Gene_Bias,
    y = prop,
    fill = dte.bias.direction
  )
) +
  geom_col(
  ) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 0.01)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  facet_wrap(~Tissue) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_fill_manual(values = c(
    "female-biased" = "#E05A4E",
    "male-biased" = "#08519c",
    "unbiased" = "#bdbdbd"
  )) +
  labs(
    x = "Differential gene expression",
    y = "Isoforms",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "plain", size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "lines")
  )

nested_props

final_plot <-
  (p_gene /
   p_isoform /
   nested_props) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right"
  )

final_plot

nested_props

################################################################################
### (6) Differential Transcript Usage (aka. differential splicing; DTU)
################################################################################

# Run fit on filtered dge
fit_g <- glmQLFit(dge_tran_gonad, 
                 design_g,
                 robust = T) # setting robust = T is strongly recommended in edgeR doc
fit_t <- glmQLFit(dge_tran_tail, 
                 design_t,
                 robust = T) 

# Using tutorial from Baldoni et al. 2025
# https://plbaldoni.rbind.io/DTU-code/
# https://github.com/plbaldoni/DTU-code

# Run diffSplice for DTU analysis
# diffSplice tests for differential usage at the transcript level grouped by genes
gonad_test_splice <- diffSplice(fit_g, 
                                contrast = contrast_matrix_g[,"Male_vs_Female"], 
                                geneid = dge_tran_gonad$genes$GeneID,      
                                exonid = dge_tran_gonad$genes$TranscriptID,
                                robust = T) 
# Total number of exons:  30844 
# Total number of genes:  10543 
# Number of genes with 1 exon:  2811 
# Mean number of exons in a gene:  3 
# Max number of exons in a gene:  50 
tail_test_splice  <- diffSplice(fit_t, 
                                contrast = contrast_matrix_t[,"Male_vs_Female"],
                                geneid = dge_tran_tail$genes$GeneID,  
                                exonid = dge_tran_tail$genes$TranscriptID,
                                robust = T) 
# Total number of exons:  25555 
# Total number of genes:  9934 
# Number of genes with 1 exon:  3215 
# Mean number of exons in a gene:  3 
# Max number of exons in a gene:  38 

# NOTE: the max exons doesn't equal the known total max (124) because
# some transcripts have 0 expression in one tissue entirely

# Each coefficient is the difference between the log-fold-change for that isoform 
# versus the average log-fold-change for all other isoforms for the same gene.
summary(gonad_test_splice$coefficients)
summary(tail_test_splice$coefficients)

# Extract results with topSplice()

# We can look at gene level and transcript level, 
# both dfs should have same genes in them
# test = "t" gives transcript-level (within-gene) tests
dtu_gonad_transcript <- topSplice(gonad_test_splice, test = "t", number = Inf)
dtu_tail_transcript <- topSplice(tail_test_splice, test = "t", number = Inf)
# test = "F" gives gene-level overall DTU tests
dtu_gonad_gene <- topSplice(gonad_test_splice, test = "F", number = Inf)
dtu_tail_gene <- topSplice(tail_test_splice, test = "F", number = Inf)

# Add significance column (FDR < 0.05 and absolute logFC >= 1)
dtu_gonad_gene$dtu.is.altspliced.gene <- ifelse(dtu_gonad_gene$FDR < 0.05, T, F)
dtu_tail_gene$dtu.is.altspliced.gene <- ifelse(dtu_tail_gene$FDR < 0.05 , T, F)

# Save
write.csv(dtu_gonad_gene,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "01_dtu_gonad_gene_long.csv"),
          row.names = F)
write.csv(dtu_tail_gene,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "02_dtu_tail_gene_long.csv"),
          row.names = F)

# Add significance column (FDR < 0.05 and absolute logFC >= 1)
dtu_gonad_transcript$dtu.is.altspliced.tran <- ifelse(dtu_gonad_transcript$FDR < 0.05, T, F)
dtu_tail_transcript$dtu.is.altspliced.tran <- ifelse(dtu_tail_transcript$FDR < 0.05 , T, F)

glimpse(dtu_gonad_gene)
table(dtu_gonad_gene$dtu.is.altspliced.gene)

glimpse(dtu_tail_gene)
table(dtu_tail_gene$dtu.is.altspliced.gene)

# Clean up
rm(dge_tran_gonad)
rm(dge_tran_tail)
rm(gonad_test_splice)
rm(tail_test_splice)
rm(contrast_matrix_g)
rm(contrast_matrix_t)
rm(design_g)
rm(design_t)

# Save DTU results
write.csv(dtu_gonad_transcript,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "01_dtu_gonad_tran_long.csv"),
          row.names = F)
write.csv(dtu_tail_transcript,
          here::here("data",
                     "isoseq",
                     "12_edgeR_dge",
                     "02_dtu_tail_tran_long.csv"),
          row.names = F)

################################################################################
### (7) Venn diagram at gene level for each tissue
################################################################################

# Create individual dfs and then join

# Select relevant columns from DGE and DTU
dge_gonad_sel <- dge_gonad_results %>%
  dplyr::select(GeneID, dge.is.biased)
dtu_gonad_sel <- dtu_gonad_gene %>%
  dplyr::select(GeneID, dtu.is.altspliced.gene)

# Full join all data by GeneID
gonad_combined_df <- full_join(dge_gonad_sel, 
                               dtu_gonad_sel, 
                               by = "GeneID")

# This will introduce NA values for dtu.is.altspliced for genes
# with only 1 transcript, let's exclude those genes for now
# Remove rows where DTU info is missing
gonad_combined_df <- gonad_combined_df %>%
  filter(!is.na(dtu.is.altspliced.gene))
gonad_combined_df$dataset <- "gonad"

# tail
dge_tail_sel <- dge_tail_results %>%
  dplyr::select(GeneID, dge.is.biased)
dtu_tail_sel <- dtu_tail_gene %>%
  dplyr::select(GeneID, dtu.is.altspliced.gene)
tail_combined_df <- full_join(dge_tail_sel, 
                              dtu_tail_sel, 
                              by = "GeneID")
tail_combined_df <- tail_combined_df %>%
  filter(!is.na(dtu.is.altspliced.gene))
tail_combined_df$dataset <- "tail"

# Combine all datasets into one dataframe
combined_df_gene <- rbind(gonad_combined_df, 
                          tail_combined_df)

# Define a list of datasets
datasets <- unique(combined_df_gene$dataset)

# Define colors for regions
region_colors <- c(
  "1"   = "#DBCB9A",  # DTU only
  "2"   = "#ECB4A0",  # DGE only
  "1/2" = "#D4AE76"   # DGE & DTU
)

# Loop through each dataset and generate Venn plot
for (dataset_name in datasets) {
  
  # Filter data for the current dataset
  current_data <- combined_df_gene %>% filter(dataset == dataset_name)
  
  # Create gene sets directly from logical columns
  gene_sets <- list(
    DTU = current_data$GeneID[current_data$dtu.is.altspliced.gene],
    DGE = current_data$GeneID[current_data$dge.is.biased]
  )
  
  # Create the Venn object
  venn_obj <- Venn(gene_sets)
  
  # Extract region data from the Venn object
  venn_data <- process_data(venn_obj)
  
  # Prepare polygon data for the sets (circles)
  set_edges <- venn_data$setEdge
  
  # Prepare polygon data for regions
  region_edges <- venn_data$regionEdge
  
  # Calculate total unique genes across all sets
  total_genes <- length(unique(current_data$GeneID))
  
  # Add a label combining count and percentage to region_labels
  region_labels <- venn_data$regionLabel %>%
    mutate(
      fill_color = region_colors[id],
      percent = round(100 * count / total_genes, 1),
      label = paste0(count, "\n(", percent, "%)")
    )
  
  # Add a label combining count and percentage to region_labels
  dge_genes <- gene_sets$DGE
  dtu_genes <- gene_sets$DTU
  biased_genes <- union(dge_genes, dtu_genes)
  all_genes <- unique(current_data$GeneID)
  unbiased_genes <- setdiff(all_genes, biased_genes)
  
  unbiased_count <- length(unbiased_genes)
  unbiased_percent <- round(100 * unbiased_count / total_genes, 1)
  unbiased_label <- paste0(unbiased_count, "\n(", unbiased_percent, "%)")
  
  # Plot
  plot_obj <- ggplot() +
    # Circle boundaries
    geom_path(data = set_edges, aes(x = X, y = Y, group = id, color = id), size = 1) +
    
    # Region fills
    geom_polygon(data = region_edges, aes(x = X, y = Y, group = id, fill = id), color = NA, alpha = 0.6) +
    
    # Region count + % labels
    geom_text(data = region_labels, aes(x = X, y = Y, label = label), size = 6, lineheight = 0.9) +
    
    # Manual fill colors
    scale_fill_manual(values = region_colors) +
    
    # Circle edge colors
    scale_color_manual(values = c("1" = "black", "2" = "black", "1/2" = "black")) +
    
    # Add label for unbiased genes (outside Venn)
    annotate("text", 
             x = min(set_edges$X) + 0.8, 
             y = max(set_edges$Y) + 0, 
             label = unbiased_label, 
             hjust = 0.5, 
             vjust = 0.5, 
             size = 6) +
    
    theme_void() +
    theme(legend.position = "none") +
    coord_fixed() +
    theme(
      plot.margin = margin(t = 30, r = 5, b = 5, l = 5)  # top, right, bottom, left in pts
    )
  
  # Save plot
  plot_filename <- paste0("08_vennplot_gene_", dataset_name, "_long.png")
  ggsave(here::here("figures", plot_filename), plot = plot_obj, width = 4, height = 6, dpi = 300)
  
  message(paste("Plot saved for", dataset_name))
}

################################################################################
### (8) Venn diagram at transcript level for each tissue
################################################################################

# Select relevant columns from DGE and DTU
dte_gonad_sel <- dte_gonad_results %>%
  dplyr::select(TranscriptID, dte.is.biased, dte.bias.direction)
dtu_gonad_sel <- dtu_gonad_transcript %>%
  dplyr::select(TranscriptID, dtu.is.altspliced.tran)

# Full join all data by GeneID
gonad_combined_df <- full_join(dte_gonad_sel, 
                               dtu_gonad_sel, 
                               by = "TranscriptID")

# This will introduce NA values for dtu.is.altspliced for genes
# with only 1 transcript, let's exclude those genes for now
# Remove rows where DTU info is missing
gonad_combined_df <- gonad_combined_df %>%
  filter(!is.na(dtu.is.altspliced.tran))
gonad_combined_df$dataset <- "gonad"

# tail
dte_tail_sel <- dte_tail_results %>%
  dplyr::select(TranscriptID, dte.is.biased, dte.bias.direction)
dtu_tail_sel <- dtu_tail_transcript %>%
  dplyr::select(TranscriptID, dtu.is.altspliced.tran)
tail_combined_df <- full_join(dte_tail_sel, 
                              dtu_tail_sel, 
                              by = "TranscriptID")
tail_combined_df <- tail_combined_df %>%
  filter(!is.na(dtu.is.altspliced.tran))
tail_combined_df$dataset <- "tail"

# Combine all datasets into one dataframe
combined_df_tran <- rbind(gonad_combined_df, tail_combined_df)

# Define a list of datasets
datasets <- unique(combined_df_tran$dataset)

# Define colors for regions
region_colors <- c(
  "1"   = "#DBCB9A",  # DTU only
  "2"   = "#CAD4C9",  # DTE only
  "1/2" = "#D2CFB2"   # DTE & DTU
)

# Loop through each dataset and generate Venn plot
for (dataset_name in datasets) {
  
  # Filter data for the current dataset
  current_data <- combined_df_tran %>% filter(dataset == dataset_name)
  
  # Create gene sets directly from logical columns
  gene_sets <- list(
    DTE = current_data$TranscriptID[current_data$dte.is.biased],
    DTU = current_data$TranscriptID[current_data$dtu.is.altspliced.tran]
  )
  
  # Create the Venn object
  venn_obj <- Venn(gene_sets)
  
  # Extract region data from the Venn object
  venn_data <- process_data(venn_obj)
  
  # Prepare polygon data for the sets (circles)
  set_edges <- venn_data$setEdge
  
  # Prepare polygon data for regions
  region_edges <- venn_data$regionEdge
  
  # Calculate total unique genes across all sets
  total_trans <- length(unique(current_data$TranscriptID))
  
  # Add a label combining count and percentage to region_labels
  region_labels <- venn_data$regionLabel %>%
    mutate(
      fill_color = region_colors[id],
      percent = round(100 * count / total_trans, 1),
      label = paste0(count, "\n(", percent, "%)")
    )
  
  # Add a label combining count and percentage to region_labels
  dte_trans <- gene_sets$DTE
  dtu_trans <- gene_sets$DTU
  biased_trans <- union(dte_trans, dtu_trans)
  all_trans <- unique(current_data$TranscriptID)
  unbiased_trans <- setdiff(all_trans, biased_trans)
  
  unbiased_count <- length(unbiased_trans)
  unbiased_percent <- round(100 * unbiased_count / total_trans, 1)
  unbiased_label <- paste0(unbiased_count, "\n(", unbiased_percent, "%)")
  
  # Plot
  plot_obj <- ggplot() +
    # Circle boundaries
    geom_path(data = set_edges, aes(x = X, y = Y, group = id, color = id), size = 1) +
    
    # Region fills
    geom_polygon(data = region_edges, aes(x = X, y = Y, group = id, fill = id), color = NA, alpha = 0.6) +
    
    # Region count + % labels
    geom_text(data = region_labels, aes(x = X, y = Y, label = label), size = 6, lineheight = 0.9) +
    
    # Manual fill colors
    scale_fill_manual(values = region_colors) +
    
    # Circle edge colors
    scale_color_manual(values = c("1" = "black", "2" = "black", "1/2" = "black")) +
    
    # Add label for unbiased genes (outside Venn)
    annotate("text", 
             x = min(set_edges$X) + 0.8, 
             y = max(set_edges$Y) + 0, 
             label = unbiased_label, 
             hjust = 0.5, 
             vjust = 0.5, 
             size = 6) +
    
    theme_void() +
    theme(legend.position = "none") +
    coord_fixed() +
    theme(
      plot.margin = margin(t = 30, r = 5, b = 5, l = 5)  # top, right, bottom, left in pts
    )
  
  # Save plot
  plot_filename <- paste0("08_vennplot_tran_", dataset_name, "_long.png")
  ggsave(here::here("figures", plot_filename), plot = plot_obj, width = 4, height = 6, dpi = 300)
  
  message(paste("Plot saved for", dataset_name))
}


