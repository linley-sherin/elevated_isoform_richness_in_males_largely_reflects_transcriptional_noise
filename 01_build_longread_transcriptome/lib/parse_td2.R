
library(Biostrings)
library(stringr)
library(dplyr)

# Load gff3 where cds are mapped back to genome ---------------------------

# Load gff file output from TD2 - these are all of the orfs that have strands matching the original genomic coordinates
gff <- read.delim(here::here("data",
                             "01_build_longread_transcriptome",
                             "13_transdecoder",
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
    gene_id     = str_extract(orf_id, "^PB\\.\\d+"),            # Extract G.10000
    transcript_id = str_extract(orf_id, "^PB\\.\\d+\\.\\d+")    # Extract G.10000.2
  )

# Load cds output and parse -----------------------------------------------

# Load cds file output from TD2
cds <- readDNAStringSet(here::here("data",
                                   "01_build_longread_transcriptome",
                                   "13_transdecoder",
                                   "transcripts.fasta.TD2.cds"),
                        format = "fasta")

# Get headers
headers <- names(cds)

# Extract different components of header
orf_id <- str_extract(headers, "^\\S+")
transcript_id <- str_extract(headers, "(?<=GENE\\.)[^~]+")
gene_id <- str_extract(transcript_id, "^[^.]+\\.[^.]+")
orf_type       <- str_extract(headers, "(?<=ORF type:)[^ ]+")
psauron_score  <- as.numeric(str_extract(headers, "(?<=psauron_score=)[0-9.]+"))
orf_len        <- as.integer(str_extract(headers, "(?<=len:)[0-9]+"))
cds_coords     <- str_extract(headers, "PB\\.[^ ]+$")

# Sequence as string
cds_seq <- as.character(cds)

# Combine into data frame
cds_df <- data.frame(
  orf_id             = orf_id,
  transcript_id      = transcript_id,
  gene_id            = gene_id,
  orf_type           = orf_type,
  psauron_score      = psauron_score,
  orf_len            = orf_len,
  cds_seq            = cds_seq,
  stringsAsFactors   = FALSE
)

# Remove all potential orfs that don't appear in gff_parsed
# This means that they do not have the correct strand data to match the genome
cds_filtered <- cds_df %>%
  filter(orf_id %in% gff_parsed$orf_id)

# Find best ORFs
cds_best <- cds_filtered %>%
  filter(psauron_score >= 0.90) %>%                     # filter by score
  filter(orf_len > 90) %>%                              # filter by length
  group_by(transcript_id) %>%                           # group by transcript
  slice_max(psauron_score, with_ties = FALSE) %>%       # keep top scoring row
  ungroup()

# Look at # unique orfs
cds_seq_counts <- cds_best %>%
  group_by(gene_id) %>%
  summarise(unique_cds_seq_count = n_distinct(cds_seq)) %>%
  arrange(desc(unique_cds_seq_count))

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

# Load ORF (AA) output and parse -----------------------------------------------

# Load AA sequences
aa <- readAAStringSet(here::here("data",
                                 "01_build_longread_transcriptome",
                                 "13_transdecoder",
                                 "transcripts.fasta.TD2.pep"),
                      format = "fasta")

# Extract headers
aa_headers <- names(aa)

# Convert to character sequences
orf_seq <- as.character(aa)

aa_df <- data.frame(
  orf_id = str_extract(aa_headers, "^\\S+"),  # the PB. ID
  orf_seq = orf_seq,
  stringsAsFactors = FALSE
)

# Merge with the rest of the TD results
cds_best_coords <- cds_best_coords %>%
  left_join(aa_df, by = "orf_id")

# Save
write.csv(cds_best_coords,
          here::here("data",
                     "01_build_longread_transcriptome",
                     "13_transdecoder",
                     "parsed_results.csv"), 
          row.names = FALSE)
