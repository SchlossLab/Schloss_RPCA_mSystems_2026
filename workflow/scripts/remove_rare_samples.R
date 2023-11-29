#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

args <-  commandArgs(trailingOnly = TRUE)
count_meta_fname <- args[1]
n_sequence_threshold <- as.numeric(args[2])

# count_meta_fname <- "data/sponges/data_treatment.tsv.gz"
# n_sequence_threshold <- 13755

count_design_fname <- str_replace(
  count_meta_fname,
  "_treatment.*",
  glue(".{n_sequence_threshold}.raw_design.tsv.gz"))

fread(count_meta_fname) %>% 
  pivot_longer(cols = starts_with("OTU_")) %>%
  mutate(n_seqs = sum(value), .by = "sample") %>%
  filter(n_seqs >= n_sequence_threshold) %>%
  select(-n_seqs) %>%
  pivot_wider(id_cols = c(sample, treatment),
              names_from = name,
              values_from = value) %>%
  write_tsv(count_design_fname)
