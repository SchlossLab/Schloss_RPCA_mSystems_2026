#!/usr/bin/env Rscript

library(tidyverse)


dir.create("data/process", showWarnings = FALSE)

list.files(path = "data/simulation", pattern = ".*knn_.*.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(.default = col_double()),
           id = "fname") %>%
  mutate(simulation = str_replace(fname,
                                  "data/simulation/(.*)\\.\\d*\\..*",
                                  "\\1"),
         depth = str_replace(fname,
                             ".*\\.(\\d*)\\..*",
                             "\\1"),
         process = str_replace(fname,
                             ".*\\.\\d*\\.(\\w*)\\.\\w*\\.knn.*",
                             "\\1"),
         calculator = str_replace(fname,
                             ".*\\.\\d*\\.\\w*\\.(\\w*)\\.knn.*",
                             "\\1"),
         approach = str_replace(fname,
                                ".*knn_(.*)\\..*",
                              "\\1"),
         process_calculator = paste(process, calculator, sep = "_") ) %>%
  select(simulation, process_calculator, approach, depth,
         iteration, k_size, test, validate) %>%
  write_tsv("data/process/simulation.knn.tsv.gz")
