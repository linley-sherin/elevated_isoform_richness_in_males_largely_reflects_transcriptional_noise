# Load edgeR results
# Note: these are already filtered so sex-limited genes are removed 

library(here)
library(tidyverse)
library(dplyr)

### Differential Gene Expression
# gonad
dge_gonad <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "01_dge_gonad_long.csv")) %>%
  dplyr::select(GeneID, logFC, dge.bias.direction) %>%
  dplyr::rename(gene = GeneID) %>%
  mutate(tissue = "gonad")
# tail
dge_tail <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "02_dge_tail_long.csv")) %>%
  dplyr::select(GeneID, logFC, dge.bias.direction) %>%
  dplyr::rename(gene = GeneID) %>%
  mutate(tissue = "tail")
# combine
combined_dge <- bind_rows(dge_gonad, dge_tail)

### Differential Transcript Expression
# gonad
dte_gonad <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "01_dte_gonad_long.csv")) %>%
  dplyr::select(GeneID, TranscriptID, logFC, dte.bias.direction) %>%
  dplyr::rename(gene = GeneID, transcript = TranscriptID) %>%
  mutate(tissue = "gonad")
# tail
dte_tail <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "02_dte_tail_long.csv")) %>%
  dplyr::select(GeneID, TranscriptID, logFC, dte.bias.direction) %>%
  dplyr::rename(gene = GeneID, transcript = TranscriptID) %>%
  mutate(tissue = "tail")
# combine
combined_dte <- bind_rows(dte_gonad, dte_tail)

### Differential Transcript Usage at Gene Level
# gonad
dtug_gonad <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "01_dtu_gonad_gene_long.csv")) %>%
  dplyr::select(GeneID, dtu.is.altspliced.gene) %>%
  dplyr::rename(gene = GeneID) %>%
  mutate(tissue = "gonad")
# tail
dtug_tail <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "02_dtu_tail_gene_long.csv")) %>%
  dplyr::select(GeneID, dtu.is.altspliced.gene) %>%
  dplyr::rename(gene = GeneID) %>%
  mutate(tissue = "tail")
# combine
combined_dtug <- bind_rows(dtug_gonad, dtug_tail)

### Differential Transcript Usage at Transcript Level
# gonad
dtut_gonad <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "01_dtu_gonad_tran_long.csv")) %>%
  dplyr::select(GeneID, TranscriptID, dtu.is.altspliced.tran) %>%
  dplyr::rename(gene = GeneID, transcript = TranscriptID) %>%
  mutate(tissue = "gonad")
# tail
dtut_tail <- read.csv(here::here("data", "isoseq", "12_edgeR_dge", "02_dtu_tail_tran_long.csv")) %>%
  dplyr::select(GeneID, TranscriptID, dtu.is.altspliced.tran) %>%
  dplyr::rename(gene = GeneID, transcript = TranscriptID) %>%
  mutate(tissue = "tail")
# combine
combined_dtut <- bind_rows(dtut_gonad, dtut_tail)
