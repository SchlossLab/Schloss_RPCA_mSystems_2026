#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/process", showWarnings = FALSE)

list.files(path = "data/simulation", pattern = ".*raw.dtest.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(method = col_character(),
                            alternative = col_character(),
                            .default = col_double()),
            id = "fnames") %>%
  mutate(simulation = str_replace(fnames,
                                  "data/simulation/(.*)\\.\\d*\\..*",
                                  "\\1"),
         depth = str_replace(fnames,
                             ".*\\.(\\d*)\\..*",
                             "\\1")) %>%
  select(simulation, depth, p_value, statistic, method, alternative) %>%
  write_tsv("data/process/simulation.dtest.tsv")
