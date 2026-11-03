# Load mastertables with metrics calculated at different levels

library(here)
library(tidyverse)
library(dplyr)

mt_samp <- read.csv(here::here("data",
                                "03_data_analysis",
                                "03_mastertable",
                                "mastertable_sample.csv"))

mt_treat <- read.csv(here::here("data",
                                "03_data_analysis",
                                "03_mastertable",
                                "mastertable_treatment.csv"))

mt_tissue <- read.csv(here::here("data",
                                "03_data_analysis",
                                "03_mastertable",
                                "mastertable_tissue.csv"))

mt_all <- read.csv(here::here("data",
                                "03_data_analysis",
                                "03_mastertable",
                                "mastertable_all.csv"))