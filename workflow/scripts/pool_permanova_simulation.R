#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/process", showWarnings = FALSE)

list.files(path = "data/simulation", pattern = ".*permanova.*.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(dmatrix_fname = col_character(),
                            .default = col_double())) %>%
  mutate(simulation = str_replace(dmatrix_fname,
                                  "data/simulation/(.*)\\.\\d*\\..*",
                                  "\\1"),
         depth = str_replace(dmatrix_fname,
                             ".*\\.(\\d*)\\..*",
                             "\\1"),
         process = str_replace(dmatrix_fname,
                             ".*\\.\\d*\\.(\\w*)\\.\\w*\\..*",
                             "\\1"),
         calculator = str_replace(dmatrix_fname,
                             ".*\\.\\d*\\.\\w*\\.(\\w*)\\..*",
                             "\\1"),
         process_calculator = paste(process, calculator, sep = "_")) %>%
  select(simulation, process_calculator, depth, f_statistic, p_value, r2) %>%
  write_tsv("data/process/simulation.permanova.tsv")
