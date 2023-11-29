#!/usr/bin/env Rscript

library(tidyverse)

dtest_fnames <- commandArgs(trailingOnly = TRUE)

# dtest_fnames <- c("data/crc/data.dtest.tsv", "data/sponges/data.dtest.tsv",
#                   "data/sleep_apnea/data.dtest.tsv")


dir.create("data/process", showWarnings = FALSE)

dtest_fnames %>%
  read_tsv(col_types = cols(method = col_character(),
                            alternative = col_character(),
                            .default = col_double()),
          id = "fnames") %>%
  mutate(case_study = str_replace(fnames,
                                  "data/(.*)/data\\.dtest\\.tsv",
                                  "\\1")) %>%
  select(case_study, depth = n_seqs, p_value, statistic, method, alternative) %>%
  write_tsv("data/process/case_studies.dtest.tsv")
