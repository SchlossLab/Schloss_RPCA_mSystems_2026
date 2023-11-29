#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/process", showWarnings = FALSE)

p <- ".*permanova.*"

fnames <- c(list.files(path = "data/crc", pattern = p, full.names = TRUE),
          list.files(path = "data/sponges", pattern = p, full.names = TRUE),
          list.files(path = "data/sleep_apnea", pattern = p, full.names = TRUE))

fnames %>%
  read_tsv(col_types = cols(dmatrix_fname = col_character(),
                            .default = col_double())) %>%
  mutate(dataset = str_replace(dmatrix_fname,
                                  "data/(.*)/data\\.\\d*\\..*",
                                  "\\1"),
         depth = str_replace(dmatrix_fname,
                             ".*\\.(\\d*)\\..*",
                             "\\1"),
         n_per_group = str_replace(dmatrix_fname,
                             ".*\\.(\\d*)_\\d*\\..*",
                             "\\1"),
         iteration = str_replace(dmatrix_fname,
                             ".*\\.\\d*_(\\d*)\\..*",
                             "\\1"),
         process = str_replace(dmatrix_fname,
                             ".*\\.\\d*\\.(\\w*)\\.\\w*\\..*",
                             "\\1"),
         calculator = str_replace(dmatrix_fname,
                             ".*\\.\\d*\\.\\w*\\.(\\w*)\\..*",
                             "\\1"),
         process_calculator = paste(process, calculator, sep = "_")) %>%
  select(dataset, n_per_group, iteration, process_calculator, depth, f_statistic, p_value, r2) %>%
  write_tsv("data/process/case_studies.permanova.tsv")
