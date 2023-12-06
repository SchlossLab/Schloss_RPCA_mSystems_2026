#!/usr/bin/env Rscript

library(tidyverse)

dtest_fnames <- commandArgs(trailingOnly = TRUE)

dir.create("data/pooled_data", showWarnings = FALSE)


dtest_fnames %>%
  read_tsv(col_types = cols(method = col_character(),
                            alternative = col_character(),
                            .default = col_double()),
          id = "fnames") %>%
  mutate(case_study = str_replace(fnames,
                                  "data/(.*)/data\\.dtest\\.tsv",
                                  "\\1")) %>%
  select(case_study, depth = n_seqs, p_value, statistic, method, alternative) %>%
  write_tsv("data/pooled_data/depth_confounding.tsv")
