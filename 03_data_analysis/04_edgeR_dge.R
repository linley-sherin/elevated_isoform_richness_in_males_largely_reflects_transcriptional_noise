## Step 04: Run edgeR for DGE analysis
## Last modified by Linley Sherin 2026-11-03
#
# Only running on genes and isoforms expressed in both sexes
#
# Runs locally.
# (1) Load files
# (2) Create design matrix
# (3) Differential Gene Expression (DGE)
# (4) Differential Transcript Expression (DTE)
# (5) Compare DGE and DTE
# (6) Differential Transcript Usage (aka. differential splicing; DTU)
# (7) Venn diagram at gene level for each tissue
# (8) Venn diagram at transcript level for each tissue
# (9) Sex-biased expression restricted to isoforms with intact ORFs
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/04_edgeR_dge.R
#
# Reads:  data/01_build_longread_transcriptome/16_edgeR_filterexp/
#           01_eR_genegonad.RData, 02_eR_genetail.RData,
#           01_eR_trangonad.RData, 02_eR_trantail.RData
#         data/03_data_analysis/03_mastertable/ (via lib/load_mt.R)
#         data/reference/metadata.csv
#         data/03_data_analysis/02_orf_status/orf_status.csv
# Writes: data/03_data_analysis/04_edgeR_dge/
#           01_dge_gonad_long.csv, 02_dge_tail_long.csv
#           01_dte_gonad_long.csv, 02_dte_tail_long.csv
#           01_dtu_gonad_gene_long.csv, 02_dtu_tail_gene_long.csv
#           01_dtu_gonad_tran_long.csv, 02_dtu_tail_tran_long.csv
#           03_dge_all_vs_coding_long.csv, 03_die_by_orfstatus_long.csv
#         figures/04_vennplot_gene_<tissue>_long.png
#         figures/04_vennplot_tran_<tissue>_long.png
#         figures/04_die_by_orfstatus.png
#         figures/04_dge_all_vs_coding_logFC.png, 04_dge_all_vs_coding_counts.png
#
# Overwrites its outputs every time it runs.

library(edgeR)
library(dplyr)
library(ggVennDiagram)
library(tidyverse)
library(limma)
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
source(here::here("03_data_analysis", "lib", "load_mt.R"))

# Get isoforms expressed in both sexes (in >=1 sample of each sex);
# sex-specific isoforms are not tested for DTE/DTU
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
                                "reference",
                                "metadata.csv"))

# Your vector of all sample names
all_samples <- c("RFG1", "RFT1", "RMG1", "RMT1",
                 "RFG2", "RFT2", "RMG2", "RMT2",
                 "RFG3", "RFT3", "RMG3", "RMT3",
                 "RFG4", "RFT4", "RMG4", "RMT4",
                 "RFG5", "RFT5", "RMG5", "RMT5",
                 "RFG6", "RFT6", "RMG6", "RMT6")

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

# Put samples in the same order as the DGE object columns (the design rows
# are matched to the count columns by position, so a different order would
# silently mislabel sex)
sample_info_gonad <- sample_info_gonad[match(colnames(dge_gene_gonad), sample_info_gonad$sample), ]
sample_info_tail <- sample_info_tail[match(colnames(dge_gene_tail), sample_info_tail$sample), ]
stopifnot(identical(sample_info_gonad$sample, colnames(dge_gene_gonad)),
          identical(sample_info_gonad$sample, colnames(dge_tran_gonad)),
          identical(sample_info_tail$sample, colnames(dge_gene_tail)),
          identical(sample_info_tail$sample, colnames(dge_tran_tail)))

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
rownames(design_t) <- sample_info_tail$sample

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
# -1*Female 1*Male
# Down               3528
# NotSig             5941
# Up                 3668
summary(decideTests(tail_test_gene))
# -1*Female 1*Male
# Down                  8
# NotSig            12471
# Up                  195

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
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "01_dge_gonad_long.csv"),
          row.names = F)
write.csv(dge_tail_results,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "02_dge_tail_long.csv"),
          row.names = F)

################################################################################
### (4) Differential Transcript Expression (DTE)
################################################################################

# Estimate dispersions with the quasi-likelihood pipeline

# Sanity check
dim(dge_tran_gonad)
dim(dge_tran_tail)

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
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "01_dte_gonad_long.csv"),
          row.names = F)
write.csv(dte_tail_results,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
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
# Total number of exons:  48055 
# Total number of genes:  13137 
# Number of genes with 1 exon:  4040 
# Mean number of exons in a gene:  4 
# Max number of exons in a gene:  88 
tail_test_splice  <- diffSplice(fit_t, 
                                contrast = contrast_matrix_t[,"Male_vs_Female"],
                                geneid = dge_tran_tail$genes$GeneID,  
                                exonid = dge_tran_tail$genes$TranscriptID,
                                robust = T) 
# Total number of exons:  44454 
# Total number of genes:  12674 
# Number of genes with 1 exon:  3971 
# Mean number of exons in a gene:  4 
# Max number of exons in a gene:  107 

# NOTE: the max exons doesn't equal the known total max (124) because
# some transcripts have 0 expression in one tissue entirely

# Each coefficient is the difference between the log-fold-change for that isoform 
# versus the average log-fold-change for all other isoforms for the same gene.

# Extract results with topSplice()

# We can look at gene level and transcript level, 
# both dfs should have same genes in them
# test = "t" gives transcript-level (within-gene) tests
dtu_gonad_transcript <- topSplice(gonad_test_splice, test = "t", number = Inf)
dtu_tail_transcript <- topSplice(tail_test_splice, test = "t", number = Inf)
# test = "F" gives gene-level overall DTU tests
dtu_gonad_gene <- topSplice(gonad_test_splice, test = "F", number = Inf)
dtu_tail_gene <- topSplice(tail_test_splice, test = "F", number = Inf)

# Add significance column (FDR < 0.05; no logFC cut-off for DTU)
dtu_gonad_gene$dtu.is.altspliced.gene <- ifelse(dtu_gonad_gene$FDR < 0.05, T, F)
dtu_tail_gene$dtu.is.altspliced.gene <- ifelse(dtu_tail_gene$FDR < 0.05 , T, F)

# Save
write.csv(dtu_gonad_gene,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "01_dtu_gonad_gene_long.csv"),
          row.names = F)
write.csv(dtu_tail_gene,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "02_dtu_tail_gene_long.csv"),
          row.names = F)

# Add significance column (FDR < 0.05; no logFC cut-off for DTU)
dtu_gonad_transcript$dtu.is.altspliced.tran <- ifelse(dtu_gonad_transcript$FDR < 0.05, T, F)
dtu_tail_transcript$dtu.is.altspliced.tran <- ifelse(dtu_tail_transcript$FDR < 0.05 , T, F)

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
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "01_dtu_gonad_tran_long.csv"),
          row.names = F)
write.csv(dtu_tail_transcript,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
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
  "1"   = "#DBCB9A",  # DTU only (set 1 = DTU)
  "2"   = "#ECB4A0",  # DGE only (set 2 = DGE)
  "1/2" = "#D4AE76"   # DGE & DTU
)

# Loop through each dataset and generate Venn plot
for (dataset_name in datasets) {
  
  # Filter data for the current dataset
  current_data <- combined_df_gene %>% filter(dataset == dataset_name)
  
  # Create gene sets directly from logical columns
  gene_sets <- list(
    DTU = current_data$GeneID[current_data$dtu.is.altspliced.gene %in% TRUE],
    DGE = current_data$GeneID[current_data$dge.is.biased %in% TRUE]
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
  plot_filename <- paste0("04_vennplot_gene_", dataset_name, "_long.png")
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
  "1"   = "#DBCB9A",  # DTE only (set 1 = DTE)
  "2"   = "#CAD4C9",  # DTU only (set 2 = DTU)
  "1/2" = "#D2CFB2"   # DTE & DTU
)

# Loop through each dataset and generate Venn plot
for (dataset_name in datasets) {
  
  # Filter data for the current dataset
  current_data <- combined_df_tran %>% filter(dataset == dataset_name)
  
  # Create gene sets directly from logical columns
  gene_sets <- list(
    DTE = current_data$TranscriptID[current_data$dte.is.biased %in% TRUE],
    DTU = current_data$TranscriptID[current_data$dtu.is.altspliced.tran %in% TRUE]
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
  plot_filename <- paste0("04_vennplot_tran_", dataset_name, "_long.png")
  ggsave(here::here("figures", plot_filename), plot = plot_obj, width = 4, height = 6, dpi = 300)
  
  message(paste("Plot saved for", dataset_name))
}

################################################################################
### (9) Sex-biased expression restricted to isoforms with intact ORFs
################################################################################

# Do isoforms without intact ORFs (premature stop/NMD,
# retained intron, non-coding) inflate estimates of sex-biased expression?
# (A) Proportion of sex-biased isoforms (DIE) by ORF status
# (B) Rerun DGE with gene counts from coding isoforms only and compare
#      with DGE from all isoforms (prevalence and magnitude of sex bias)

# Intact ORF = ptc_status "coding" (complete ORF, no NMD signal, no retained intron)
orf_status <- read.csv(here::here("data",
                                  "03_data_analysis",
                                  "02_orf_status",
                                  "orf_status.csv"))[,-1]

coding_tx <- orf_status$transcript[orf_status$ptc_status == "coding"]

# Colours used throughout
bias_colors <- c(
  "female-biased" = "#E05A4E",
  "unbiased"      = "#bdbdbd",
  "male-biased"   = "#08519c"
)

### (A) Proportion of sex-biased isoforms by ORF status

die_orf <- bind_rows(
  dte_gonad_results %>% mutate(Tissue = "Gonad"),
  dte_tail_results %>% mutate(Tissue = "Tail")
) %>%
  left_join(dplyr::select(orf_status, transcript, ptc_status),
            by = c("TranscriptID" = "transcript")) %>%
  mutate(
    ptc_status = recode(ptc_status,
                        "coding" = "Coding",
                        "nmd" = "Premature STOP",
                        "intron_retention" = "Retained intron",
                        "non-coding" = "Non-coding"),
    ptc_status = factor(ptc_status, levels = c("Coding", "Premature STOP",
                                               "Retained intron", "Non-coding")),
    Tissue = factor(Tissue, levels = c("Tail", "Gonad")),
    dte.bias.direction = factor(dte.bias.direction, levels = names(bias_colors))
  )

# Sanity check: every tested isoform should have an ORF status
cat("Isoforms without ORF status:", sum(is.na(die_orf$ptc_status)), "\n")

die_orf_summary <- die_orf %>%
  count(Tissue, ptc_status, dte.bias.direction, .drop = FALSE) %>%
  group_by(Tissue, ptc_status) %>%
  mutate(n_isoforms = sum(n),
         prop = n / n_isoforms) %>%
  ungroup()

print(die_orf_summary %>% filter(dte.bias.direction != "unbiased"), n = Inf)

# Are isoforms without intact ORFs more often sex-biased than coding ones?
for (t in levels(die_orf$Tissue)) {
  x <- die_orf %>% filter(Tissue == t)
  cat("\n---", t, ": sex-biased (DIE) vs intact ORF ---\n")
  print(chisq.test(table(intact_orf = x$ptc_status == "Coding",
                         biased = x$dte.is.biased)))
}

# Plot: % female-biased, unbiased and male-biased isoforms per ORF status.
# Number of isoforms tested is shown under each category.
die_orf_plot <- die_orf_summary %>%
  mutate(category = paste0(sub(" ", "\n", ptc_status), "\n(n = ", scales::comma(n_isoforms), ")"),
         category = factor(category, levels = unique(category[order(ptc_status)])),
         label = ifelse(prop < 0.01,
                        scales::percent(prop, accuracy = 0.01),
                        scales::percent(prop, accuracy = 0.1)))

p_die_orf <- ggplot(die_orf_plot,
                    aes(x = category, y = prop, fill = dte.bias.direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.6) +
  facet_wrap(~Tissue, scales = "free") +
  scale_fill_manual(values = bias_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "DIE results (% of isoforms)") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey90", color = NA))

p_die_orf

ggsave(here::here("figures", "04_die_by_orfstatus.png"),
       plot = p_die_orf, width = 9, height = 4.5, dpi = 300)

### (B) DGE with coding isoforms only vs. all isoforms

# Reload transcript-level DGE objects (unfiltered; removed after DTU above).
# These contain every isoform of the genes kept in step 16, so summing them
# gives the same gene counts used for the original DGE.
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "01_eR_trangonad.RData"))
load(here::here("data",
                "01_build_longread_transcriptome",
                "16_edgeR_filterexp",
                "02_eR_trantail.RData"))

# Same expression filter as step 16: TPM > 1 in >= 3 samples of a sex
# and total TPM > 1
filter_by_treatment <- function(dge_obj, TPM_cutoff = 1, min_samples = 3, total_cutoff = 1) {
  TPM_matrix <- cpm(dge_obj)
  keep_by_group <- rep(FALSE, nrow(TPM_matrix))
  for (grp in unique(dge_obj$samples$treatment)) {
    idx <- which(dge_obj$samples$treatment == grp)
    keep_by_group <- keep_by_group | rowSums(TPM_matrix[, idx, drop = FALSE] > TPM_cutoff) >= min_samples
  }
  keep_by_group & rowSums(TPM_matrix) > total_cutoff
}

# Sum isoform counts to genes (only isoforms in keep_tx), normalise, and run
# the same DGE as section (3) on the genes in genes_tested
run_dge <- function(dge_tran, keep_tx, genes_tested, apply_filter) {
  keep <- rownames(dge_tran) %in% keep_tx
  gene_counts <- rowsum(dge_tran$counts[keep, , drop = FALSE],
                        group = dge_tran$genes$GeneID[keep])
  dge <- DGEList(counts = gene_counts,
                 samples = dge_tran$samples[, c("sample", "sex", "tissue", "treatment")])
  dge$genes <- data.frame(GeneID = rownames(dge), row.names = rownames(dge))
  dge <- dge[rowSums(dge$counts) > 0, , keep.lib.sizes = FALSE]
  if (apply_filter) dge <- dge[filter_by_treatment(dge), , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  # Test only the genes from the original DGE (lib sizes kept from all genes, as in (1))
  dge <- dge[rownames(dge) %in% genes_tested, ]

  design <- model.matrix(~ 0 + sex, data = dge$samples)
  colnames(design) <- c("Female", "Male")
  contrast <- makeContrasts(Male_vs_Female = Male - Female, levels = colnames(design))

  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  test <- glmTreat(fit, contrast = contrast[, "Male_vs_Female"], lfc = 1)
  res <- edgeR::topTags(test, n = Inf)$table

  # Same bias calls as section (3)
  res$dge.is.biased <- res$FDR < 0.05 & abs(res$logFC) >= 1
  res$dge.bias.direction <- "unbiased"
  res$dge.bias.direction[res$dge.is.biased & res$logFC > 1] <- "male-biased"
  res$dge.bias.direction[res$dge.is.biased & res$logFC < -1] <- "female-biased"
  res
}

# Check: all isoforms should reproduce the original DGE
dge_gonad_check <- run_dge(dge_tran_gonad, rownames(dge_tran_gonad), dge_gonad_results$GeneID, apply_filter = FALSE)
dge_tail_check <- run_dge(dge_tran_tail, rownames(dge_tran_tail), dge_tail_results$GeneID, apply_filter = FALSE)
cat("Gonad, all isoforms rebuilt:\n"); print(table(dge_gonad_check$dge.bias.direction))
cat("Tail, all isoforms rebuilt:\n"); print(table(dge_tail_check$dge.bias.direction))

# Coding isoforms only (genes need >= 1 coding isoform and must pass the
# expression filter on coding counts alone)
dge_gonad_coding <- run_dge(dge_tran_gonad, coding_tx, dge_gonad_results$GeneID, apply_filter = TRUE)
dge_tail_coding <- run_dge(dge_tran_tail, coding_tx, dge_tail_results$GeneID, apply_filter = TRUE)

# Share of each gene's counts that comes from coding isoforms, per sex
prop_coding <- function(dge_tran) {
  is_coding <- rownames(dge_tran) %in% coding_tx
  by_sex <- lapply(split(seq_len(ncol(dge_tran)), dge_tran$samples$sex), function(idx) {
    all_counts <- rowsum(rowSums(dge_tran$counts[, idx, drop = FALSE]), dge_tran$genes$GeneID)
    coding_counts <- rowsum(rowSums(dge_tran$counts[is_coding, idx, drop = FALSE]),
                            dge_tran$genes$GeneID[is_coding])
    out <- setNames(rep(0, nrow(all_counts)), rownames(all_counts))
    out[rownames(coding_counts)] <- coding_counts[, 1]
    out / all_counts[, 1]
  })
  data.frame(GeneID = names(by_sex$Female),
             prop_coding_female = by_sex$Female,
             prop_coding_male = by_sex$Male[names(by_sex$Female)])
}

# Join all-isoform and coding-only results
compare_dge <- function(all_res, coding_res, dge_tran, tissue) {
  all_res %>%
    dplyr::select(GeneID, logFC_all = logFC, FDR_all = FDR, bias_all = dge.bias.direction) %>%
    left_join(coding_res %>%
                dplyr::select(GeneID, logFC_coding = logFC, FDR_coding = FDR,
                              bias_coding = dge.bias.direction),
              by = "GeneID") %>%
    left_join(prop_coding(dge_tran), by = "GeneID") %>%
    mutate(Tissue = tissue,
           tested_coding = !is.na(bias_coding))
}

dge_compare <- bind_rows(
  compare_dge(dge_gonad_results, dge_gonad_coding, dge_tran_gonad, "Gonad"),
  compare_dge(dge_tail_results, dge_tail_coding, dge_tran_tail, "Tail")
) %>%
  mutate(Tissue = factor(Tissue, levels = c("Tail", "Gonad")),
         bias_all = factor(bias_all, levels = names(bias_colors)),
         bias_coding = factor(bias_coding, levels = names(bias_colors)))

# Genes that could not be tested with coding isoforms only
cat("\nGenes not tested with coding isoforms only (no coding isoform or below filter):\n")
print(dge_compare %>% count(Tissue, bias_all, tested_coding) %>% filter(!tested_coding))

# Prevalence: % sex-biased genes, on genes tested both ways
dge_shared <- dge_compare %>% filter(tested_coding)

prevalence <- dge_shared %>%
  pivot_longer(c(bias_all, bias_coding), names_to = "analysis", values_to = "bias") %>%
  mutate(analysis = recode(analysis, "bias_all" = "All isoforms",
                           "bias_coding" = "Coding isoforms only")) %>%
  count(Tissue, analysis, bias, .drop = FALSE) %>%
  group_by(Tissue, analysis) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

print(prevalence, n = Inf)

# Which genes change category (rows = all isoforms, columns = coding only)
for (t in levels(dge_shared$Tissue)) {
  x <- dge_shared %>% filter(Tissue == t)
  cat("\n---", t, ": DGE all isoforms (rows) vs coding only (columns) ---\n")
  tab <- table(all = x$bias_all, coding = x$bias_coding)
  print(tab)
  print(round(prop.table(tab, margin = 1) * 100, 2))
}

# Magnitude: |logFC| all vs coding only
magnitude <- dge_shared %>%
  group_by(Tissue) %>%
  summarise(
    n_genes = n(),
    r_logFC = cor(logFC_all, logFC_coding),
    median_abs_logFC_all = median(abs(logFC_all)),
    median_abs_logFC_coding = median(abs(logFC_coding)),
    wilcox_p = wilcox.test(abs(logFC_coding), abs(logFC_all), paired = TRUE)$p.value,
    .groups = "drop"
  )
print(magnitude)

# Same for genes sex-biased with all isoforms
magnitude_biased <- dge_shared %>%
  filter(bias_all != "unbiased") %>%
  group_by(Tissue, bias_all) %>%
  summarise(
    n_genes = n(),
    median_abs_logFC_all = median(abs(logFC_all)),
    median_abs_logFC_coding = median(abs(logFC_coding)),
    median_prop_coding_female = median(prop_coding_female),
    median_prop_coding_male = median(prop_coding_male),
    .groups = "drop"
  )
print(magnitude_biased)

# Save
write.csv(dge_compare,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "03_dge_all_vs_coding_long.csv"),
          row.names = F)
write.csv(die_orf_summary,
          here::here("data",
                     "03_data_analysis",
                     "04_edgeR_dge",
                     "03_die_by_orfstatus_long.csv"),
          row.names = F)

# Plot: logFC with all isoforms vs coding isoforms only
logfc_data <- dge_shared %>%
  mutate(change = case_when(
    bias_all == bias_coding ~ "Same",
    bias_all == "unbiased" ~ "Biased only with coding",
    bias_coding == "unbiased" ~ "Biased only with all isoforms",
    TRUE ~ "Changed direction"
  ),
  change = factor(change, levels = c("Same", "Biased only with all isoforms",
                                     "Biased only with coding", "Changed direction"))) %>%
  arrange(change)  # draw "Same" first so the changed genes sit on top

change_colors <- c("Same" = "grey75",
                   "Biased only with all isoforms" = "#E69F00",
                   "Biased only with coding" = "#009E73",
                   "Changed direction" = "black")

p_logfc <- ggplot(logfc_data, aes(x = logFC_all, y = logFC_coding, color = change)) +
  geom_point(size = 0.9, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-1, 1), color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), color = "grey60", linewidth = 0.3) +
  facet_wrap(~Tissue) +
  scale_color_manual(values = change_colors, name = NULL, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Male:female log2FC (all isoforms)",
       y = "Male:female log2FC (coding isoforms only)") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey90", color = NA))

p_logfc

ggsave(here::here("figures", "04_dge_all_vs_coding_logFC.png"),
       plot = p_logfc, width = 8, height = 4, dpi = 300)

# Plot: number of sex-biased genes, all vs coding only
p_prev <- prevalence %>%
  filter(bias != "unbiased") %>%
  ggplot(aes(x = analysis, y = n, fill = bias)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = scales::comma(n)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3) +
  facet_wrap(~Tissue, scales = "free_y") +
  scale_fill_manual(values = bias_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Sex-biased genes") +
  theme_classic() +
  theme(strip.background = element_rect(fill = "grey90", color = NA))

p_prev

ggsave(here::here("figures", "04_dge_all_vs_coding_counts.png"),
       plot = p_prev, width = 7, height = 4, dpi = 300)
