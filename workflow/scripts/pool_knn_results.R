#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)

knn_results_fnames <- commandArgs(trailingOnly = TRUE)

knn_results_fnames %>%
  read_tsv(col_types = cols(.default = col_double()),
           id = "fname") %>%
  separate(fname, sep = "[\\/\\.]",
          into = c("data_dir", "dataset", "data", "depth",
                  "size_rep", "process", "distance", "knn_approach", "tsv")) %>%
  separate(size_rep, into = c("nsamples_per_group", "rep")) %>%
  select(dataset, depth, nsamples_per_group, rep, process, distance, knn_approach,
         iteration, k_size, test, validate) %>%
  write_tsv("data/pooled_data/knn_results.tsv.gz")
