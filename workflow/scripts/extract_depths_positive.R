#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

simulation_file_name <- args[1]
depth <- as.numeric(args[2])

# depth <- 2000

simulation_file_name %>%
  read_csv(col_types = cols(OTUs = col_character(),
                            .default = col_double())) %>%
  filter(sequence_depth == depth) %>%
  select(OTUs, starts_with("sample")) %>%
  column_to_rownames("OTUs") %>%
  t() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(glue("data/simulation/positive_raw.{depth}.tsv.gz"))
