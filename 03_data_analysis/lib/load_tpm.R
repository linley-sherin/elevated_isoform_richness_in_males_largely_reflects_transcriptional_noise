# Load tpm matrices with metrics calculated at different levels

library(tidyverse)
library(dplyr)
library(here)

# Load TPM matrices for gonad and tail
tpm_gonad <- read.csv(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "16_edgeR_filterexp",
                                 "01_gonadtran_TMMTPMmatrix.csv"))
tpm_tail <- read.csv(here::here("data",
                                "01_build_longread_transcriptome",
                                "16_edgeR_filterexp",
                                "02_tailtran_TMMTPMmatrix.csv"))[-1]

# Reshape TPM matrices and apply them to each tissue
reshape_tpm <- function(df, tissue_name) {
  df %>%
    pivot_longer(
      cols = matches("^X\\d"),
      names_to = "sample",
      values_to = "isoformTPM"
    ) %>%
    mutate(
      sex = case_when(
        str_detect(sample, "F") ~ "female",
        str_detect(sample, "M") ~ "male",
        TRUE ~ NA_character_
      ),
      tissue = tissue_name,
      log2_isoformTPM = log2(isoformTPM + 1)) %>%  # log2 transform per sample
    group_by(sample, gene) %>%
    mutate(
      geneTPM = sum(isoformTPM),
      log2_geneTPM = (log2(geneTPM + 1)),
      ra_sample = isoformTPM/geneTPM, # relative abundance within this sample
      ra_sample = if_else(is.nan(ra_sample), 0, ra_sample) # Set NaN values to 0
    ) %>%
    ungroup()
}

# Apply to counts and bind rows
gonad_long <- reshape_tpm(tpm_gonad, "gonad")
tail_long <- reshape_tpm(tpm_tail, "tail")
combined_samp <- bind_rows(gonad_long, tail_long)
