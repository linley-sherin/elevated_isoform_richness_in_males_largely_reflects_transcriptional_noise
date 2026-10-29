### This script looks at isoforms that have different directions of sex bias
### to get at biologically meaningful splice variants between the sexes.

# (1) Load
# (2) Find putative isoform switches
# (3) Categorize switch type
# (4) Classify ORF outcomes independent of type
# (5) Summary table for plotting
# (6) Plot

library(ggplot2)
library(patchwork)
library(ltc)

################################################################################
### (1) Load
################################################################################

source("load_dge.R")
source("load_mt.R")

# Clean up files you don't need
rm(dge_gonad,
   dge_tail,
   dte_gonad,
   dte_tail,
   dtug_gonad,
   dtug_tail,
   combined_dtug,
   dtut_gonad,
   dtut_tail,
   combined_dtut,
   mt_samp,
   mt_tissue,
   mt_all)

# Join mt with edgeR results and orf_status column
mt_treat <- mt_treat %>%
  left_join(combined_dte, by = c("gene", "transcript", "tissue")) %>%
  left_join(combined_dge, by = c("gene", "tissue"))

################################################################################
### (2) Find putative isoform switches
################################################################################

# Find isoform switching genes
isoform_switching_genes <- mt_treat %>%
  # Keep only sex-biased isoforms
  filter(dte.bias.direction != "unbiased") %>%
  # Keep only effective isoforms
  filter(effGene_treat == TRUE) %>%
  # CRITICAL UPDATE: Ensure >5 TPM is in the sex matching the bias direction
  filter(
    (dte.bias.direction == "male-biased" & sex == "male" & mean_isoformTPM_treat > 5) |
      (dte.bias.direction == "female-biased" & sex == "female" & mean_isoformTPM_treat > 5)
  ) %>%
  # Check that gene has transcripts biased in opposite directions
  group_by(gene, tissue) %>%
  filter("male-biased" %in% dte.bias.direction && "female-biased" %in% dte.bias.direction) %>%
  ungroup() %>%
  dplyr::select(gene, transcript, tissue, sex, mean_isoformTPM_treat, mean_geneTPM_treat, ptc_status, orf_identity, dte.bias.direction)

# Count the total number of unique switching genes per tissue
tissue_counts <- isoform_switching_genes %>% 
  group_by(tissue) %>% 
  summarise(unique_genes = n_distinct(gene), .groups = "drop")

# Output the counts for each tissue using cat
purrr::walk2(tissue_counts$tissue, tissue_counts$unique_genes, function(tish, count) {
  cat("Total number of unique isoform switching genes in", tish, ":", count, "\n")
})

# Count the number of transcripts per gene
transcripts_per_gene <- isoform_switching_genes %>%
  group_by(gene) %>%
  summarize(num_transcripts = n_distinct(transcript)) %>%
  arrange(desc(num_transcripts)) # Sort with the highest transcript counts first

# View the resulting dataframe
print(transcripts_per_gene)

#### Do isoform switching genes also produce unbiased isoforms?

# Look for whether isoform switching genes also have unbiased isoforms expressed >5TPM in both sexes
unbiased_isoforms <- mt_treat %>%
  filter(
    paste(gene, tissue) %in%
      paste(isoform_switching_genes$gene, isoform_switching_genes$tissue),
    dte.bias.direction == "unbiased",
    mean_isoformTPM_treat > 5
  )

# Number of unique isoform-switching genes
n_total <- isoform_switching_genes %>%
  distinct(gene, tissue) %>%
  nrow()

# Number with ≥1 unbiased isoform >5 TPM
n_with <- unbiased_isoforms %>%
  distinct(gene, tissue) %>%
  nrow()

# Number without
n_without <- n_total - n_with

id_summary <- tibble(
  category = c("With \nunbiased isoforms", "Without \nunbiased isoforms"),
  n_genes = c(n_with, n_without)
) %>%
  mutate(
    category = factor(
      category,
      levels = c("Without \nunbiased isoforms", "With \nunbiased isoforms")
    ),
    proportion = n_genes / sum(n_genes),
    label = n_genes
  )

ggplot(id_summary, aes(x = category, y = proportion)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = label), vjust = -0.5) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1)
  ) +
  labs(
    x = NULL,
    y = "Proportion of \nisoform-switching genes"
  ) +
  theme_classic()

### How many unbiased isoforms have orfs?
unbiased_orf_summary <- unbiased_isoforms %>%
  distinct(gene, tissue, transcript, orf_status, orf_identity) %>%
  summarise(
    n_isoforms = n(),
    n_with_orf = sum(!is.na(orf_identity)),
    n_without_orf = sum(is.na(orf_identity)),
    prop_with_orf = n_with_orf / n_isoforms
  )

unbiased_orf_summary

### What is there orf_status
unbiased_isoforms %>%
  distinct(gene, tissue, transcript, orf_status) %>%
  count(orf_status)

### Do they code the same ORF as the sex-biased isoforms
switching_orfs <- isoform_switching_genes %>%
  filter(!is.na(orf_identity)) %>%
  distinct(gene, tissue, orf_identity) %>%
  mutate(shared_with_switch = TRUE)

comparison <- unbiased_isoforms %>%
  distinct(gene, tissue, transcript, orf_identity) %>%
  left_join(
    switching_orfs,
    by = c("gene", "tissue", "orf_identity")
  ) %>%
  mutate(
    shared_with_switch = tidyr::replace_na(shared_with_switch, FALSE)
  )

comparison %>%
  filter(!is.na(orf_identity)) %>%
  count(shared_with_switch)

comparison %>%
  group_by(gene, tissue) %>%
  summarise(
    has_unbiased_unique_orf =
      any(!shared_with_switch & !is.na(orf_identity)),
    .groups = "drop"
  ) %>%
  count(has_unbiased_unique_orf)

orf_summary <- comparison %>%
  group_by(gene, tissue) %>%
  summarise(
    has_orf = any(!is.na(orf_identity)),
    has_unique_orf = any(!shared_with_switch & !is.na(orf_identity)),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      !has_orf ~ "No ORF",
      has_unique_orf ~ "Unique ORF",
      TRUE ~ "Shared ORF"
    )
  )

n_no_orf <- sum(orf_summary$category == "No ORF")
n_shared <- sum(orf_summary$category == "Shared ORF")
n_unique <- sum(orf_summary$category == "Unique ORF")

plot_df <- tibble(
  x = c(
    "Without\nunbiased isoforms",
    rep("With\nunbiased isoforms", 3)
  ),
  fill = c(
    "No unbiased isoforms",
    "No ORF",
    "Shared ORF",
    "Unique ORF"
  ),
  n = c(
    n_without,
    n_no_orf,
    n_shared,
    n_unique
  )
) %>%
  mutate(proportion = n / n_total)

ggplot(plot_df, aes(x = x, y = proportion, fill = fill)) +
  geom_col(width = 0.9) +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5), 
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.8)
  ) +
  scale_fill_manual(
    values = c(
      "No unbiased isoforms" = "grey70",
      "No ORF" = "#83AF9B",
      "Shared ORF" = "#f8da8a",
      "Unique ORF" = "#8bc8cb"
    ),
    breaks = c("No ORF", "Shared ORF", "Unique ORF")
  ) +
  labs(
    x = NULL,
    y = "Proportion of\nisoform-switching genes",
    fill = "Unbiased isoforms"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

################################################################################
### (3) Categorize switch type
################################################################################

switch_summary <- isoform_switching_genes %>%
  group_by(tissue, gene) %>%
  summarise(
    n_male_biased = sum(dte.bias.direction == "male-biased"),
    n_female_biased = sum(dte.bias.direction == "female-biased"),
    .groups = "drop"
  ) %>%
  mutate(
    switch_type = case_when(
      n_male_biased == 1 & n_female_biased == 1 ~ "1Female:1Male",
      n_male_biased > 1 & n_female_biased == 1  ~ "1Female:MultiMale",
      n_male_biased == 1 & n_female_biased > 1  ~ "MultiFemale:1Male",
      TRUE ~ "MultiFemale:MultiMale"
    )
  )

################################################################################
### (4) Classify ORF outcomes independent of type
################################################################################

# Add switch architecture back onto isoform-level data
isoform_switching_genes2 <- isoform_switching_genes %>%
  left_join(
    switch_summary %>% 
      dplyr::select(tissue, gene, switch_type),
    by = c("tissue", "gene")
  )

orf_switch_summary <- isoform_switching_genes2 %>%
  group_by(tissue, gene, switch_type) %>%
  summarise(
    male_orfs = list(unique(
      orf_identity[dte.bias.direction == "male-biased"]
    )),
    female_orfs = list(unique(
      orf_identity[dte.bias.direction == "female-biased"]
    )),
    .groups = "drop"
  ) %>%
  mutate(
    # remove non-coding isoforms for coding comparisons
    male_coding = map(male_orfs, ~ .x[!is.na(.x)]),
    female_coding = map(female_orfs, ~ .x[!is.na(.x)]),
    n_male_coding = map_int(male_coding, length),
    n_female_coding = map_int(female_coding, length),
    shared = map2(
      male_coding,
      female_coding,
      intersect
    ),
    n_shared = map_int(shared, length),
    male_unique = map2(
      male_coding,
      shared,
      setdiff
    ),
    female_unique = map2(
      female_coding,
      shared,
      setdiff
    ),
    n_male_unique = map_int(male_unique, length),
    n_female_unique = map_int(female_unique, length),
    male_has_nc = map_lgl(
      male_orfs,
      ~ any(is.na(.x))
    ),
    female_has_nc = map_lgl(
      female_orfs,
      ~ any(is.na(.x))
    )
  ) %>%
  mutate(
    orf_outcome = case_when(
      # no coding ORFs in either sex
      n_male_coding == 0 &
        n_female_coding == 0 ~
        "All non-coding",
      # coding lost in one sex
      n_male_coding == 0 |
        n_female_coding == 0 ~
        "Loss of coding in one sex",
      # identical coding repertoires and no non-coding isoforms
      n_shared > 0 &
        n_male_unique == 0 &
        n_female_unique == 0 &
        !male_has_nc &
        !female_has_nc ~
        "All redundant ORFs",
      # completely different coding repertoires and no non-coding isoforms
      n_shared == 0 &
        !male_has_nc &
        !female_has_nc ~
        "All sex-biased ORFs",
      # all other combinations
      TRUE ~
        "Mix"
    )
  )

# Check when loss of coding occurs, what sex it's in 
orf_switch_summary %>%
  filter(orf_outcome == "Loss of coding in one sex") %>%
  mutate(
    loss_direction = case_when(
      n_male_coding == 0 & n_female_coding > 0 ~ "Male loss \nof coding",
      n_female_coding == 0 & n_male_coding > 0 ~ "Female loss \nof coding"
    )
  ) %>%
  count(tissue, loss_direction)

loss_df <- orf_switch_summary %>%
  filter(orf_outcome == "Loss of coding in one sex") %>%
  mutate(
    loss_direction = case_when(
      n_male_coding == 0 & n_female_coding > 0 ~ "Male loss \nof coding",
      n_female_coding == 0 & n_male_coding > 0 ~ "Female loss \nof coding"
    )
  )

loss_prop <- loss_df %>%
  count(loss_direction) %>%
  mutate(
    proportion = n / sum(n),
    label = n
  )

ggplot(loss_prop, aes(x = loss_direction, y = proportion)) +
  geom_col(fill = "#C88E82") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5), 
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.1)
  ) +
  labs(
    x = NULL,
    y = "Proportion of \nisoform-switching genes"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    )
  )

# counts
loss_counts <- table(loss_df$loss_direction)

loss_counts

binom.test(
  x = loss_counts["Male loss \nof coding"],
  n = sum(loss_counts),
  p = 0.5
)

################################################################################
### (5) Summary table for plotting
################################################################################

orf_plot_df <- orf_switch_summary %>%
  count(
    tissue,
    switch_type,
    orf_outcome,
    name = "Count"
  ) %>%
  group_by(tissue, switch_type) %>%
  mutate(
    Percentage = round(
      Count / sum(Count) * 100,
      1
    ),
    label = paste0(Percentage, "%")
  ) %>%
  ungroup()

################################################################################
### (6) Plot
################################################################################

orf_plot_df <- orf_plot_df %>%
  mutate(
    orf_outcome = factor(
      orf_outcome,
      levels = c(
        "All sex-biased ORFs",
        "All redundant ORFs",
        "Loss of coding in one sex",
        "All non-coding",
        "Mix"
      )
    )
  )

coastal_palette <- c(
  "All sex-biased ORFs"        = "#8A9A86",
  "All redundant ORFs"         = "#7FA1B4",
  "Loss of coding in one sex"  = "#C88E82",
  "All non-coding"            = "#D1B48C",
  "Mix"                        = "#A394A8"
)

final_plot <- ggplot(
  orf_plot_df,
  aes(
    y = "",
    x = Count,
    fill = orf_outcome
  )
) +
  geom_col(
    width = 1,
    #color = "white",
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(
    values = coastal_palette,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = scales::comma_format(),
    expand = c(0,0)
  ) +
  facet_wrap(
    ~switch_type,
    ncol = 1
  ) +
  labs(
    x = "Isoform-switching genes",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(
      fill = "#F8F4EF",
      color = NA
    ),
    panel.background = element_rect(
      fill = "#F8F4EF",
      color = NA
    ),
    legend.background = element_rect(
      fill = "#F8F4EF",
      color = NA
    ),
    legend.box.background = element_rect(
      fill = "#F8F4EF",
      color = NA
    ),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0, "lines"),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(final_plot)
