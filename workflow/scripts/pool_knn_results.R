#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)

n_fnames <- as.numeric(commandArgs(trailingOnly = TRUE))

knn_results_fnames <- list.files(path = "data",
                                 pattern = ".*knn_.*.tsv",
                                 recursive = TRUE,
                                 full.names = TRUE)

stopifnot(length(knn_results_fnames) == n_fnames)

knn_results_fnames %>%
  read_tsv(col_types = cols(.default = col_double()),
           id = "fname") %>%
  separate(fname, sep = "[\\/\\.]",
          into = c("data_dir", "dataset", "data", "depth",
                  "size_rep", "process", "distance", "knn_approach", "tsv")) %>%
  separate(size_rep, into = c("nsamples_per_group", "rep")) %>%
  select(dataset, depth, nsamples_per_group, rep, process, distance,
         knn_approach, iteration, k_size, test, validate) %>%
  write_tsv("data/pooled_data/knn_results.tsv.gz")
