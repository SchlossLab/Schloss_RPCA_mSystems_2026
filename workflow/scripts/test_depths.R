#!/usr/bin/env Rscript

library(tidyverse)
library(broom)

args <-  commandArgs(trailingOnly = TRUE)

tsv_fname <- args[1]
design_fname <- args[2]
output_file_name <- args[3]

# tsv_fname <- "data/simulation/overlap_raw.10000.tsv.gz"
# design_fname <- "data/simulation/overlap_raw.design.gz"

read_tsv(tsv_fname, col_types = cols(sample = col_character(),
                                     .default = col_double())) %>%
  select(sample, starts_with("OTU")) %>%
  pivot_longer(cols = starts_with("OTU")) %>%
  group_by(sample) %>%
  summarize(n = sum(value), .groups = "drop") %>%
  inner_join(., read_tsv(design_fname, col_types = "cc"), by = "sample") %>%
  wilcox.test(n ~ treatment, data = .) %>%
  tidy() %>%
  rename(p_value = p.value) %>%
  write_tsv(output_file_name)
