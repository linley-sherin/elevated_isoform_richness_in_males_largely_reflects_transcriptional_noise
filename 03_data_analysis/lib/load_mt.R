# Load mastertables with metrics calculated at different levels

library(here)
library(tidyverse)
library(dplyr)

mt_samp <- read.csv(here::here("data",
                                "isoseq",
                                "11_mastertable",
                                "mastertable_sample.csv"))

mt_treat <- read.csv(here::here("data",
                                "isoseq",
                                "11_mastertable",
                                "mastertable_treatment.csv"))

mt_tissue <- read.csv(here::here("data",
                                "isoseq",
                                "11_mastertable",
                                "mastertable_tissue.csv"))

mt_all <- read.csv(here::here("data",
                                "isoseq",
                                "11_mastertable",
                                "mastertable_all.csv"))