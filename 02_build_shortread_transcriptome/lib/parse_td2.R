
library(Biostrings)
library(stringr)
library(dplyr)

# Load gff3 where cds are mapped back to genome ---------------------------

# Load gff file output from TD2 - these are all of the orfs that have strands matching the original genomic coordinates
gff <- read.delim(here::here("data",
                             "02_build_shortread_transcriptome",
                             "02_transdecoder",
                             "transcripts_with_cds.gff3"), 
                  header = FALSE, 
                  comment.char = "#", 
                  stringsAsFactors = FALSE)
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Extract relevant attributes
gff_parsed <- gff %>%
  filter(type == "CDS") %>%
  mutate(
    parent      = str_extract(attributes, "Parent=([^;]+)"),   # Extract Parent
    orf_id      = str_replace(parent, "Parent=", ""),          # This is the real ORF ID
    transcript = str_extract(orf_id, ".*(?=\\.p\\d+$)")
  )

# Load t2g map
t2g <- rtracklayer::import(here::here("data",
                                      "02_build_shortread_transcriptome",
                                      "01_stringtie",
                                      "merge",
                                      "merged.gtf")) %>%
  data.frame() %>%
  dplyr::select(gene_id, ref_gene_id, transcript_id) %>%
  dplyr::rename(transcript = transcript_id, gene = gene_id, ref_gene = ref_gene_id) %>%
  filter(!is.na(gene)) %>%
  distinct()

glimpse(t2g)

# Parse the GFF attributes and join
gff_cleaned <- gff %>%
  mutate(
    transcript = str_extract(attributes, 'Name=([^ ]+)') %>% 
      str_replace('Name=', '')
  ) %>%
  left_join(t2g, by = "transcript") %>%
  filter(type == "gene") %>%
  mutate(
    # Change "gene" to "transcript" in the type column
    type = if_else(type == "gene", "transcript", type)) %>%
  distinct()

glimpse(gff_cleaned)

# Load cds output and parse -----------------------------------------------

# Load cds file output from TD2
cds <- readDNAStringSet(here::here("data",
                                   "02_build_shortread_transcriptome",
                                   "02_transdecoder",
                                   "transcripts.fasta.TD2.cds"),
                        format = "fasta")

# Extract headers and sequences
headers <- names(cds)
orf_seq <- as.character(cds)

# Create tibble and parse metadata
cds_df <- tibble(
  header = headers,
  orf_seq = orf_seq
) %>%
  mutate(
    orf_id = str_extract(header, "^\\S+"),
    transcript = str_remove(orf_id, "\\.p\\d+$"),  # clean transcript IDs
    orf_type = str_extract(header, "(?<=ORF type:)\\S+"),
    psauron_score = as.numeric(str_extract(header, "(?<=psauron_score=)[0-9.]+")),
    orf_len = as.integer(str_extract(header, "(?<=len:)\\d+")),
    orf_coords = str_extract(header, "(?<=:)[0-9]+-[0-9]+")
  )

glimpse(cds_df)

# Remove all potential orfs that don't appear in gff_parsed
# This means that they do not have the correct strand data to match the genome
cds_filtered <- cds_df %>%
  filter(orf_id %in% gff_parsed$orf_id)

# Find best ORFs
cds_best <- cds_filtered %>%
  filter(psauron_score >= 0.90) %>%                     # Step 1: filter by score
  filter(orf_len > 90) %>%                              # Step 2: filter by length
  group_by(transcript) %>%                           # Step 3: group by transcript
  slice_max(psauron_score, with_ties = FALSE) %>%       # keep top scoring row
  ungroup()

glimpse(cds_best)

# Join with t2g map using reference gene
cds_best <- cds_best %>%
  left_join(t2g %>% select(transcript, ref_gene), by = "transcript") %>%
  rename(gene = ref_gene) %>%
  filter(!is.na(gene))

glimpse(cds_best)

# Look at # unique orfs
orf_seq_counts <- cds_best %>%
  group_by(gene) %>%
  summarise(unique_orf_seq_count = n_distinct(orf_seq)) %>%
  arrange(desc(unique_orf_seq_count))

# Look at distribution of orf_type
table(cds_best$orf_type)

# Find genomic coordinates of cds region
cds_coords <- gff_parsed %>%
  filter(type == "CDS") %>%
  group_by(orf_id, seqid) %>%
  summarise(
    cds_start = min(start),
    cds_end   = max(end),
    .groups = "drop"
  )

# Join
cds_best_coords <- cds_best %>%
  left_join(cds_coords, by = "orf_id")

# Save
write.csv(cds_best_coords,
          here::here("data",
                     "02_build_shortread_transcriptome",
                     "02_transdecoder",
                     "parsed_results.csv"), 
          row.names = FALSE)


