#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

simulation_csv_file <- args[1]
output_file <- args[2]

sample_ids <- simulation_csv_file %>%
  read_csv(col_types = cols(OTUs = col_character(),
                            .default = col_double())) %>%
  select(starts_with("sample")) %>%
  colnames()

tibble(sample = sample_ids,
       treatment = rep(c("A", "B"), each = 100)) %>%
  write_tsv(output_file)
