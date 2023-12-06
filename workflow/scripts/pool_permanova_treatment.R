#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)


parse_fnames <- function(x){
 # x <- "data/sponges/data.1000.all_1.gemelli_3.euclidean.permanova.tsv"

  read_tsv(x,
          col_types = cols(.default = col_character())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data_dir", "dataset", "data", "depth", "samples_rep",
                      "process", "distance", "dist", "gz")) %>%
    separate(samples_rep, into = c("nsamples_per_group", "rep")) %>%
    select(dataset, depth, nsamples_per_group, rep, process, distance,
          f_statistic, p_value, r2)

}

permanova_fnames <- commandArgs(trailingOnly = TRUE)

permanova_fnames %>%
  map_dfr(., parse_fnames) %>%
  write_tsv("data/pooled_data/permanova_treatment.tsv")
