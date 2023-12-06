#!/usr/bin/env Rscript

library(tidyverse)


dir.create("data/process", showWarnings = FALSE)

list.files(path = c("data/sponges", "data/crc", "data/sleep_apnea"),
    pattern = ".*knn_.*.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(.default = col_double()),
           id = "fname") %>%
  separate(fname, sep = "[\\/\\.]",
          into = c("data_dir", "dataset", "data", "depth",
                  "size_rep", "process", "distance", "knn_approach", "tsv")) %>%
  separate(size_rep, into = c("nsamples", "rep")) %>%
  select(dataset, depth, nsamples, rep, process, distance, knn_approach,
         iteration, k_size, test, validate) %>%
  write_tsv("data/process/case_studies.knn.tsv.gz")
