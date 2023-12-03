#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)


args <-  commandArgs(trailingOnly = TRUE)
count_meta_fname <- args[1]
seed <- args[2]

set.seed(seed)

# count_meta_fname <- "data/sponges/data.13755.raw_design.tsv.gz"
# seed <- 1

count_fname <- str_replace(count_meta_fname,
                          "\\.\\d*\\.raw_design.tsv.gz",
                          glue(".5000.size_{seed}.raw.tsv.gz"))

design_fname <- str_replace(count_meta_fname,
                          "\\.\\d*\\.raw_design.tsv.gz",
                            glue(".5000.size_{seed}.design"))

count_design <- fread(count_meta_fname) %>%
  pivot_longer(cols = starts_with("OTU")) %>%
  filter(value > 0) %>%
  uncount(weights = value)

count_design_10k <- count_design %>%
  slice_sample(n = 10000, by = c(sample, treatment), replace = FALSE) %>%
  count(sample, treatment, name) %>%
  mutate(sample = paste0(sample, "_10k"))

count_design_5k <- count_design %>%
  slice_sample(n = 5000, by = c(sample, treatment), replace = FALSE) %>%
  count(sample, treatment, name) %>%
  mutate(sample = paste0(sample, "_5k"))

cound_design_pool <- bind_rows(count_design_10k, count_design_5k) %>%
  pivot_wider(names_from = name, values_from = n, values_fill = 0)


cound_design_pool %>%
  select(sample, starts_with("OTU")) %>%
  write_tsv(count_fname)

cound_design_pool %>%
  select(sample, treatment) %>%
  write_tsv(design_fname)
