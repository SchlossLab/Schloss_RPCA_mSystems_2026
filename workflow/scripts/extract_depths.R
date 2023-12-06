#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
count_fname <- args[1]
output_fname <- args[2]

# count_fname <- "data/sponges/data_treatment.tsv.gz"
# output_fname <- "data/sponges/data.depths.tsv"

fread(count_fname,
      colClasses = c(sample = "character")) %>%
  pivot_longer(cols = -c(sample, treatment)) %>%
  group_by(sample, treatment) %>%
  summarize(n_seqs = sum(value)) %>%
  write_tsv(output_fname)
