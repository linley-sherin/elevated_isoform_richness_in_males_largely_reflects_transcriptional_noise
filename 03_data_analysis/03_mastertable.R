### Create master table by calculating isoform metrics at different levels

# (1) Load tpm matrices + orf_status
# (2) Calculate mean tpm, EV (expression variance) within conditions + relAbund
# (3) Calculate isoform richness (potential), diversity (perplexity)
# (4) Save master table(s)

################################################################################
### (1) Load tpm matrices + orf_status
################################################################################

source("load_tpm.R")

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


# Load sqanti output and orf status
sqanti_orf <- read.csv(here::here("data",
                                  "01_build_longread_transcriptome",
                                  "16_edgeR_filterexp",
                                  "sqanti_w_orf.csv"))[,-1]

################################################################################
### (2) Calculate mean tpm, EV (expression variance) within conditions + relAbund
################################################################################

# Summarize Gene metrics independently
gene_metrics <- combined_samp %>%
  # Get one row per gene per sample to calculate true means
  distinct(sample, gene, sex, tissue, geneTPM, log2_geneTPM) %>%
  group_by(gene) %>%
  mutate(
    mean_geneTPM_all      = mean(geneTPM, na.rm = TRUE),
    log2_mean_geneTPM_all = log2(mean_geneTPM_all + 1),
    cv_geneTPM_all        = sd(log2_geneTPM, na.rm = TRUE) / mean(log2_geneTPM, na.rm = TRUE)
  ) %>%
  group_by(gene, tissue) %>%
  mutate(
    mean_geneTPM_tissue      = mean(geneTPM, na.rm = TRUE),
    log2_mean_geneTPM_tissue = log2(mean_geneTPM_tissue + 1),
    cv_geneTPM_tissue        = sd(log2_geneTPM, na.rm = TRUE) / mean(log2_geneTPM, na.rm = TRUE)
  ) %>%
  group_by(gene, tissue, sex) %>%
  mutate(
    mean_geneTPM_treat      = mean(geneTPM, na.rm = TRUE),
    log2_mean_geneTPM_treat = log2(mean_geneTPM_treat + 1),
    cv_geneTPM_treat        = sd(log2_geneTPM, na.rm = TRUE) / mean(log2_geneTPM, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Keep only the columns we need to join back
  dplyr::select(sample, gene, contains("mean_geneTPM"), contains("cv_geneTPM"))

# Calculate Isoform metrics and merge with Gene metrics
combined_samp <- combined_samp %>%
  # Drop the old gene-level mean columns if they exist to avoid .x / .y conflicts
  dplyr::select(-contains("mean_geneTPM"), -contains("cv_geneTPM")) %>%
  
  group_by(transcript) %>%
  mutate(
    mean_isoformTPM_all      = mean(isoformTPM, na.rm = TRUE),
    log2_mean_isoformTPM_all = log2(mean_isoformTPM_all + 1),
    cv_isoformTPM_all        = sd(log2_isoformTPM, na.rm = TRUE) / mean(log2_isoformTPM, na.rm = TRUE)
  ) %>%
  group_by(transcript, tissue) %>%
  mutate(
    mean_isoformTPM_tissue      = mean(isoformTPM, na.rm = TRUE),
    log2_mean_isoformTPM_tissue = log2(mean_isoformTPM_tissue + 1),
    cv_isoformTPM_tissue        = sd(log2_isoformTPM, na.rm = TRUE) / mean(log2_isoformTPM, na.rm = TRUE)
  ) %>%
  group_by(transcript, tissue, sex) %>%
  mutate(
    mean_isoformTPM_treat      = mean(isoformTPM, na.rm = TRUE),
    log2_mean_isoformTPM_treat = log2(mean_isoformTPM_treat + 1),
    cv_isoformTPM_treat        = sd(log2_isoformTPM, na.rm = TRUE) / mean(log2_isoformTPM, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Join the robust gene metrics back
  left_join(gene_metrics, by = c("sample", "gene")) %>%
  # 3. Calculate RA levels safely
  mutate(
    ra_all    = if_else(mean_geneTPM_all > 0, mean_isoformTPM_all / mean_geneTPM_all, 0),
    ra_tissue = if_else(mean_geneTPM_tissue > 0, mean_isoformTPM_tissue / mean_geneTPM_tissue, 0),
    ra_treat  = if_else(mean_geneTPM_treat > 0, mean_isoformTPM_treat / mean_geneTPM_treat, 0)
  )

# Join with orf_status column of sqanti_w_orf
master_table <- combined_samp %>%
  left_join(dplyr::select(sqanti_orf, orf_status, gene, transcript, ptc_status, cds_identity, orf_identity), by = c("gene", "transcript"))

################################################################################
### (3) Calculate isoform richness (potential), diversity (perplexity)
################################################################################

# For this part, split master_table into sample level, treatment level
# and tissue level to make it go faster, then join at end

## Sample level
mt_samp <- master_table %>%
  dplyr::select(gene,
                transcript,
                sample,
                tissue,
                sex,
                isoformTPM,
                log2_isoformTPM,
                geneTPM,
                log2_geneTPM,
                ra_sample,
                orf_status,
                ptc_status,
                cds_identity,
                orf_identity)

# Calculate perplexity at the sample level
mt_samp <- mt_samp %>%
  filter(geneTPM > 0, isoformTPM > 0) %>%
  group_by(sample, gene) %>%
  mutate(
    # Calculate potential, entropy and perplexity at the gene level
    potential_samp = n_distinct(transcript, na.rm = TRUE),
    hGene_sample = -sum(ra_sample * log2(ra_sample), na.rm = TRUE),
    pGene_sample = 2^hGene_sample,
    # # pGene_sample = if_else(geneTPM == 0, 0, pGene_sample),
    # # Calculate AS rate
    # major = isoformTPM == max(isoformTPM),  # no need for ifelse here
    # as_sample = sum(ra_sample[!major], na.rm = TRUE),
    # # Calculate tau (isoform specificity)
    # tau_sample = {
    #   n <- n()              # scalar
    #   if (n > 1) {
    #     sum(1 - isoformTPM / max(isoformTPM, na.rm = TRUE), na.rm = TRUE) / (n - 1)
    #     } else {
    #       NA_real_
    #     }
    #   },
    # Calculate the number of effective isoforms based on pGene (rounded)
    n_eff_gene = round(pGene_sample),
    # Rank isoforms by ra_sample and label top n_eff_gene as effective
    ### NOTE: FOR ALL LEVELS I BREAK TIES WITH 'FIRST' EVEN THOUGH THIS MAY ARBITRAILTY PICK AN ISOFORM AT LOW TPM
    effGene_sample = ifelse(rank(-ra_sample, ties.method = "first") <= n_eff_gene, TRUE, FALSE)) %>%
  dplyr::select(-hGene_sample)

### Calculate protein coding level metrics
ptc_metrics <- mt_samp %>%
  dplyr::filter(ptc_status == "coding") %>%
  group_by(sample, gene) %>%
  mutate(
    # Calculate potential, entropy and perplexity at the PTC level
    ra_PTC = isoformTPM / sum(isoformTPM, na.rm = TRUE),
    hPTC_sample = -sum(ra_PTC * log2(ra_PTC), na.rm = TRUE),
    pPTC_sample = 2^hPTC_sample) %>%
  dplyr::select(-hPTC_sample)

# Join the calculated PTC metrics back to mt_samp + mark effective
mt_samp <- mt_samp %>%
  left_join(
    dplyr::select(
      ptc_metrics,
      c("sample", "gene", "transcript", "ra_PTC", "pPTC_sample")
    ),
    by = c("sample", "gene", "transcript")
  ) %>%
  group_by(sample, gene) %>%
  mutate(
    ra_PTC = if_else(is.na(ra_PTC) | is.nan(ra_PTC), 0, ra_PTC),
    # SUB IN NA VALUES WITH PERPLEXITY FROM OTHER TRANSCRIPTS
    pPTC_sample = {
      x <- pPTC_sample
      if (all(is.na(x) | is.nan(x))) 0
      else dplyr::first(x[!is.na(x) & !is.nan(x)])
    }, 
    # Calculate the number of effective isoforms based on pPTC (rounded)
    n_eff_ptc = round(pPTC_sample),
    # Rank isoforms by ra_PTC and label top n_eff_ptc as effective
    effPTC_sample = rank(-ra_PTC, ties.method = "first") <= n_eff_ptc
  ) %>%
  ungroup()

### Calculate orf-level metrics
orf_metrics <- mt_samp %>%
  # dplyr::filter(effGene_sample == TRUE) %>% # Keep only transcripts effective at the gene level
  dplyr::select(sample, gene, cds_identity, orf_identity, isoformTPM, ptc_status) %>%
  group_by(sample, gene, orf_identity) %>%
  mutate(
    # Sum isoform TPM by ORF identity (coding isoforms only)
    orf_isoformTPM = sum(isoformTPM[ptc_status == "coding"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Collapse transcripts into ORF-level rows
  distinct(sample, gene, orf_identity, .keep_all = TRUE) %>%
  dplyr::select(sample, gene, orf_identity, orf_isoformTPM) %>%
  group_by(sample, gene) %>%
  mutate(
    # Relative abundance per ORF
    ra_ORF = orf_isoformTPM / sum(orf_isoformTPM, na.rm = TRUE),
    ra_ORF = if_else(is.nan(ra_ORF), 0, ra_ORF), # Set NaN values to 0
    # ORF-level entropy and perplexity
    hORF_sample = {
      if (sum(ra_ORF, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_ORF * log2(ra_ORF), na.rm = TRUE)
      }
    },
    pORF_sample = {
      if (sum(ra_ORF, na.rm = TRUE) == 0) {
        0
      } else {
        2^hORF_sample[1]
      }
    },
    # Calculate the number of effective isoforms based on pORF (rounded)
    n_eff_orf = round(pORF_sample),
    # Rank isoforms by ra_ORF and label top n_eff_orf as effective
    efforf_sample = ifelse(rank(-ra_ORF, ties.method = "first") <= n_eff_orf, TRUE, FALSE)
  ) %>%
  ungroup() %>%
 dplyr::select(sample, gene, orf_identity, ra_ORF, pORF_sample, n_eff_orf, efforf_sample)

# Join the calculated ORF metrics back to mt_samp + mark effective
mt_samp <- mt_samp %>%
  left_join(orf_metrics, by = c("sample", "gene", "orf_identity"))

## Treatment level
mt_treat <- master_table %>%
  dplyr::select(gene,
                transcript,
                tissue,
                sex,
                mean_isoformTPM_treat,
                log2_mean_isoformTPM_treat,
                cv_geneTPM_treat,
                mean_geneTPM_treat,
                log2_mean_geneTPM_treat,
                cv_isoformTPM_treat,
                ra_treat,
                orf_status,
                ptc_status,
                cds_identity,
                orf_identity) %>%
  distinct()

# Calculate metrics at the treat level
mt_treat <- mt_treat %>%
  group_by(sex, tissue, gene) %>%
  filter(mean_geneTPM_treat > 0, mean_isoformTPM_treat > 0) %>%
  mutate(
    # Calculate potential, entropy and perplexity at the gene level
    potential_treat = n_distinct(transcript[mean_isoformTPM_treat > 0], na.rm = TRUE),   
    hGene_treat = -sum(ra_treat * log2(ra_treat), na.rm = TRUE),
    pGene_treat = 2^hGene_treat,
    # Calculate AS rate
    major = mean_isoformTPM_treat == max(mean_isoformTPM_treat),  # no need for ifelse here
    as_treat = sum(ra_treat[!major], na.rm = TRUE),
    # Calculate tau (isoform specificity)
    tau_treat = {
      n <- n()
      if (n > 1) {
        sum(1 - mean_isoformTPM_treat[mean_isoformTPM_treat > 0] / max(mean_isoformTPM_treat, na.rm = TRUE), na.rm = TRUE) / (n - 1)
      } else {
        NA_real_
      }
    },
    # Calculate potential, entropy and perplexity at the PTC level
    ra_PTC_treat = ifelse(ptc_status == "coding",
                    mean_isoformTPM_treat / sum(mean_isoformTPM_treat[ptc_status == "coding"], na.rm = TRUE),
                    0),
    hPTC_treat = {
      if (sum(ra_PTC_treat, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_PTC_treat * log2(ra_PTC_treat), na.rm = TRUE)
      }
    },
    pPTC_treat = {
      if (sum(ra_PTC_treat, na.rm = TRUE) == 0) {
        0
      } else {
        2^hPTC_treat[1]
      }
    }
  ) %>%
  ungroup() %>%
  dplyr::select(-hGene_treat, -hPTC_treat)

# Classify isoforms as effective based on pGene + pPTC
mt_treat <- mt_treat %>%
  group_by(gene, sex, tissue) %>%  # Group by gene, sex, and tissue
  mutate(
    # Calculate the number of effective isoforms based on pGene (rounded)
    n_eff_gene = round(pGene_treat),
    # Rank isoforms by ra_treat within each gene-sex-tissue group, and label the top n_eff_gene as effective
    effGene_treat = ifelse(rank(-ra_treat, ties.method = "first") <= n_eff_gene, TRUE, FALSE),
    # Calculate the number of effective isoforms based on pPTC (rounded)
    n_eff_ptc = round(pPTC_treat),
    # Rank isoforms by ra_PTC within each gene-sex-tissue group
    effPTC_treat = ifelse(rank(-ra_PTC_treat, ties.method = "first") <= n_eff_ptc, TRUE, FALSE)
  ) %>%
  ungroup()  # Ungroup after the operation is done

# Calculate ORF-level metrics
orf_metrics_treat <- mt_treat %>%
  group_by(sex, tissue, gene) %>%
  
  dplyr::select(sex, tissue, gene, orf_identity, mean_isoformTPM_treat, ptc_status) %>%
  group_by(sex, tissue, gene, orf_identity) %>%
  mutate(
    # Sum isoform TPM by orf_identity for ORF-level calculations
    orf_mean_isoformTPM_treat = sum(mean_isoformTPM_treat[ptc_status == "coding"], na.rm = TRUE)  # Sum isoform TPM for coding isoforms only
  ) %>%
  ungroup() %>%
  distinct(sex, tissue, gene, orf_identity, .keep_all = TRUE) %>%  # Remove duplicates by collapsing transcripts into ORF-level
  dplyr::select(sex, tissue, gene, orf_identity, orf_mean_isoformTPM_treat) %>%  # Keep relevant columns
  group_by(sex, tissue, gene) %>%
  mutate(
    # Calculate relative abundance per ORF identity
    ra_ORF_treat = orf_mean_isoformTPM_treat / sum(orf_mean_isoformTPM_treat, na.rm = TRUE),
    ra_ORF_treat = if_else(is.nan(ra_ORF_treat), 0, ra_ORF_treat), # Set NaN values to 0
    # Calculate entropy (hORF_treat) and perplexity (pORF_treat) for ORF level
    hORF_treat = {
      if (sum(ra_ORF_treat, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_ORF_treat * log2(ra_ORF_treat), na.rm = TRUE)
      }
    },
    pORF_treat = {
      if (sum(ra_ORF_treat, na.rm = TRUE) == 0) {
        0
      } else {
        2^hORF_treat[1]
      }
    },
    n_eff_orf = round(pORF_treat),
    efforf_treat = rank(-ra_ORF_treat, ties.method = "first") <= n_eff_orf) %>%
  ungroup() %>%
  dplyr::select(sex, tissue, gene, orf_identity, ra_ORF_treat, pORF_treat, n_eff_orf, efforf_treat)  # Keep relevant columns for ORF metrics

# Join ORF-level metrics back into mt_treat
mt_treat <- mt_treat %>%
  left_join(orf_metrics_treat, by = c("sex", "tissue", "gene", "orf_identity"))

## Tissue level
mt_tissue <- master_table %>%
  dplyr::select(gene,
                transcript,
                tissue,
                mean_isoformTPM_tissue,
                log2_mean_isoformTPM_tissue,
                cv_geneTPM_tissue,
                mean_geneTPM_tissue,
                log2_mean_geneTPM_tissue,
                cv_isoformTPM_tissue,
                ra_tissue,
                orf_status,
                ptc_status,
                cds_identity,
                orf_identity) %>%
  distinct()

# Calculate metrics at the tissue level
mt_tissue <- mt_tissue %>%
  group_by(tissue, gene) %>%
  filter(mean_geneTPM_tissue > 0, mean_isoformTPM_tissue > 0) %>%
  mutate(
    # Calculate potential, entropy and perplexity at the gene level
    potential_tissue = n_distinct(transcript[mean_isoformTPM_tissue > 0], na.rm = TRUE),   
    hGene_tissue = -sum(ra_tissue * log2(ra_tissue), na.rm = TRUE),
    pGene_tissue = 2^hGene_tissue,
    # Calculate AS rate
    major = mean_isoformTPM_tissue == max(mean_isoformTPM_tissue),  # no need for ifelse here
    as_tissue = sum(ra_tissue[!major], na.rm = TRUE),
    # Calculate tau (isoform specificity)
    tau_tissue = {
      n <- n()              # scalar
      if (n > 1) {
        sum(1 - mean_isoformTPM_tissue[mean_isoformTPM_tissue > 0] / max(mean_isoformTPM_tissue, na.rm = TRUE), na.rm = TRUE) / (n - 1)
      } else {
        NA_real_
      }
    },
    # Calculate potential, entropy and perplexity at the PTC level
    ra_PTC_tissue = ifelse(ptc_status == "coding",
                          mean_isoformTPM_tissue / sum(mean_isoformTPM_tissue[ptc_status == "coding"], na.rm = TRUE),
                          0),
    hPTC_tissue = {
      if (sum(ra_PTC_tissue, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_PTC_tissue * log2(ra_PTC_tissue), na.rm = TRUE)
      }
    },
    pPTC_tissue = {
      if (sum(ra_PTC_tissue, na.rm = TRUE) == 0) {
        0
      } else {
        2^hPTC_tissue[1]
      }
    }
  ) %>%
  ungroup() %>%
  # Calculate potential, entropy and perplexity at the orf_level
  group_by(tissue, gene, orf_identity) %>%
  mutate(
    orf_mean_isoformTPM_tissue = ifelse(ptc_status == "coding",
                                       sum(mean_isoformTPM_tissue[ptc_status == "coding"], na.rm = TRUE),
                                       NA),
    orf_mean_isoformTPM_tissue = ifelse(is.na(orf_identity), NA, orf_mean_isoformTPM_tissue)
  ) %>%
  ungroup() %>%
  dplyr::select(-hGene_tissue, -hPTC_tissue)

# Classify isoforms as effective based on pGene, pPTC, and pORF
mt_tissue <- mt_tissue %>%
  group_by(gene, tissue) %>%  # Group by gene and tissue
  mutate(
    # Calculate the number of effective isoforms based on pGene (rounded)
    n_eff_gene = round(pGene_tissue),
    # Rank isoforms by ra_tissue and label top n_eff_gene as effective
    effGene_tissue = ifelse(rank(-ra_tissue, ties.method = "first") <= n_eff_gene, TRUE, FALSE),
    # Calculate the number of effective isoforms based on pPTC (rounded)
    n_eff_ptc = round(pPTC_tissue),
    # Rank isoforms by ra_PTC and label top n_eff_ptc as effective
    effPTC_tissue = ifelse(rank(-ra_PTC_tissue, ties.method = "first") <= n_eff_ptc, TRUE, FALSE)
  ) %>%
  ungroup()

# Calculate ORF-level metrics
orf_metrics_tissue <- mt_tissue %>%
  dplyr::select(tissue, gene, orf_identity, mean_isoformTPM_tissue, ptc_status) %>%
  group_by(tissue, gene, orf_identity) %>%
  mutate(
    # Sum isoform TPM by orf_identity for ORF-level calculations
    orf_mean_isoformTPM_tissue = sum(mean_isoformTPM_tissue[ptc_status == "coding"], na.rm = TRUE)  # Sum isoform TPM for coding isoforms only
  ) %>%
  ungroup() %>%
  distinct(tissue, gene, orf_identity, .keep_all = TRUE) %>%  # Remove duplicates by collapsing transcripts into ORF-level
  dplyr::select(tissue, gene, orf_identity, orf_mean_isoformTPM_tissue) %>%  # Keep relevant columns
  group_by(tissue, gene) %>%
  mutate(
    # Calculate relative abundance per ORF identity
    ra_ORF_tissue = orf_mean_isoformTPM_tissue / sum(orf_mean_isoformTPM_tissue, na.rm = TRUE),
    # Calculate entropy (hORF_tissue) and perplexity (pORF_tissue) for ORF level
    hORF_tissue = {
      if (sum(ra_ORF_tissue, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_ORF_tissue * log2(ra_ORF_tissue), na.rm = TRUE)
      }
    },
    pORF_tissue = {
      if (sum(ra_ORF_tissue, na.rm = TRUE) == 0) {
        0
      } else {
        2^hORF_tissue[1]
      }
    },
    # Calculate the number of effective isoforms based on pORF (rounded)
    n_eff_orf = round(pORF_tissue),
    # Rank isoforms by ra_ORF and label top n_eff_orf as effective
    efforf_tissue = ifelse(rank(-ra_ORF_tissue, ties.method = "first") <= n_eff_orf, TRUE, FALSE)
  ) %>%
  ungroup() %>%
  dplyr::select(tissue, gene, orf_identity, ra_ORF_tissue, pORF_tissue, n_eff_orf, efforf_tissue)  # Keep relevant columns for ORF metrics

# Join ORF-level metrics back into mt_tissue
mt_tissue <- mt_tissue %>%
  left_join(orf_metrics_tissue, by = c("tissue", "gene", "orf_identity"))

## All data
mt_all <- master_table %>%
  dplyr::select(gene,
                transcript,
                mean_isoformTPM_all,
                log2_mean_isoformTPM_all,
                cv_geneTPM_all,
                mean_geneTPM_all,
                log2_mean_geneTPM_all,
                cv_isoformTPM_all,
                ra_all,
                orf_status,
                ptc_status,
                cds_identity,
                orf_identity) %>%
  distinct()

# Calculate metrics
mt_all <- mt_all %>%
  group_by(gene) %>%
  filter(mean_geneTPM_all > 0, mean_isoformTPM_all > 0) %>%
  mutate(
    # Calculate potential, entropy and perplexity at the gene level
    potential_all = n_distinct(transcript[mean_isoformTPM_all > 0], na.rm = TRUE),   
    hGene_all = -sum(ra_all * log2(ra_all), na.rm = TRUE),
    pGene_all = 2^hGene_all,
    # Calculate AS rate
    major = mean_isoformTPM_all == max(mean_isoformTPM_all),  # no need for ifelse here
    as_all = sum(ra_all[!major], na.rm = TRUE),
    # Calculate tau (isoform specificity)
    tau_all = {
      n <- n()              # scalar
      if (n > 1) {
        sum(1 - mean_isoformTPM_all[mean_isoformTPM_all > 0] / max(mean_isoformTPM_all, na.rm = TRUE), na.rm = TRUE) / (n - 1)
      } else {
        NA_real_
      }
    },
    # Calculate potential, entropy and perplexity at the PTC level
    ra_PTC_all = ifelse(ptc_status == "coding",
                           mean_isoformTPM_all / sum(mean_isoformTPM_all[ptc_status == "coding"], na.rm = TRUE),
                           0),
    hPTC_all = {
      if (sum(ra_PTC_all, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_PTC_all * log2(ra_PTC_all), na.rm = TRUE)
      }
    },
    pPTC_all = {
      if (sum(ra_PTC_all, na.rm = TRUE) == 0) {
        0
      } else {
        2^hPTC_all[1]
      }
    }
  ) %>%
  ungroup() %>%
  # Calculate potential, entropy and perplexity at the orf_level
  group_by(gene, orf_identity) %>%
  mutate(
    orf_mean_isoformTPM_all = ifelse(ptc_status == "coding",
                                        sum(mean_isoformTPM_all[ptc_status == "coding"], na.rm = TRUE),
                                        NA),
    orf_mean_isoformTPM_all = ifelse(is.na(orf_identity), NA, orf_mean_isoformTPM_all)
  ) %>%
  ungroup() %>%
  dplyr::select(-hGene_all, -hPTC_all)

# Classify isoforms as effective based on pGene, pPTC, and pORF
mt_all <- mt_all %>%
  group_by(gene) %>%  # Group by gene
  mutate(
    # Calculate the number of effective isoforms based on pGene (rounded)
    n_eff_gene = round(pGene_all),
    # Rank isoforms by ra_all and label top n_eff_gene as effective
    effGene_all = ifelse(rank(-ra_all, ties.method = "first") <= n_eff_gene, TRUE, FALSE),
    # Calculate the number of effective isoforms based on pPTC (rounded)
    n_eff_ptc = round(pPTC_all),
    # Rank isoforms by ra_PTC and label top n_eff_ptc as effective
    effPTC_all = ifelse(rank(-ra_PTC_all, ties.method = "first") <= n_eff_ptc, TRUE, FALSE)
  ) %>%
  ungroup()

# Calculate ORF-level metrics
orf_metrics_all <- mt_all %>%
  dplyr::select(gene, orf_identity, mean_isoformTPM_all, ptc_status) %>%
  group_by(gene, orf_identity) %>%
  mutate(
    # Sum isoform TPM by orf_identity for ORF-level calculations
    orf_mean_isoformTPM_all = sum(mean_isoformTPM_all[ptc_status == "coding"], na.rm = TRUE)  # Sum isoform TPM for coding isoforms only
  ) %>%
  ungroup() %>%
  distinct(gene, orf_identity, .keep_all = TRUE) %>%  # Remove duplicates by collapsing transcripts into ORF-level
  dplyr::select(gene, orf_identity, orf_mean_isoformTPM_all) %>%  # Keep relevant columns
  group_by(gene) %>%
  mutate(
    # Calculate relative abundance per ORF identity
    ra_ORF_all = orf_mean_isoformTPM_all / sum(orf_mean_isoformTPM_all, na.rm = TRUE),
    # Calculate entropy (hORF_all) and perplexity (pORF_all) for ORF level
    hORF_all = {
      if (sum(ra_ORF_all, na.rm = TRUE) == 0) {
        0
      } else {
        -sum(ra_ORF_all * log2(ra_ORF_all), na.rm = TRUE)
      }
    },
    pORF_all = {
      if (sum(ra_ORF_all, na.rm = TRUE) == 0) {
        0
      } else {
        2^hORF_all[1]
      }
    },
    # Calculate the number of effective isoforms based on pORF (rounded)
    n_eff_orf = round(pORF_all),
    # Rank isoforms by ra_ORF and label top n_eff_orf as effective
    efforf_all = ifelse(rank(-ra_ORF_all, ties.method = "first") <= n_eff_orf, TRUE, FALSE)
  ) %>%
  ungroup() %>%
  dplyr::select(gene, orf_identity, ra_ORF_all, pORF_all, n_eff_orf, efforf_all)  # Keep relevant columns for ORF metrics

# Join ORF-level metrics back into mt_all
mt_all <- mt_all %>%
  left_join(orf_metrics_all, by = c("gene", "orf_identity"))

################################################################################
### (4) Save master table(s)
################################################################################

write.csv(mt_samp,
          here::here("data",
                     "isoseq",
                     "11_perplexity",
                     "mastertable_sample.csv"),
          row.names = F)

write.csv(mt_treat,
          here::here("data",
                     "isoseq",
                     "11_perplexity",
                     "mastertable_treatment.csv"),
          row.names = F)

write.csv(mt_tissue,
          here::here("data",
                     "isoseq",
                     "11_perplexity",
                     "mastertable_tissue.csv"),
          row.names = F)

write.csv(mt_all,
          here::here("data",
                     "isoseq",
                     "11_perplexity",
                     "mastertable_all.csv"),
          row.names = F)

