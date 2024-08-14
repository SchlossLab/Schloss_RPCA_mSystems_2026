#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)

depth_fnames <- commandArgs(trailingOnly = TRUE)

read_tsv(depth_fnames, id = "dataset") %>%
  mutate(dataset = str_replace(dataset, "data/(.*)/data.depths.tsv", "\\1")) %>%
  write_tsv("data/pooled_data/case_study_ns.tsv")
 