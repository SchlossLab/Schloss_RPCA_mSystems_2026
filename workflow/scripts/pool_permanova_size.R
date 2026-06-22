#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)


parse_fnames <- function(x){
 # x <- "data/sponges/data.5000.size_1.gemelli_3.euclidean.permanova.tsv"

  read_tsv(x,
          col_types = cols(.default = col_character())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data_dir", "dataset", "data", "depth", "samples_rep",
                      "process", "distance", "dist", "gz")) %>%
    separate(samples_rep, into = c("samples", "rep")) %>%
    select(dataset, depth, rep, process, distance, effect,
          f_statistic, p_value, r2)

}

permanova_fnames <- commandArgs(trailingOnly = TRUE)

permanova_fnames %>%
  map_dfr(., parse_fnames) %>%
  write_tsv("data/pooled_data/permanova_size.tsv.gz")
