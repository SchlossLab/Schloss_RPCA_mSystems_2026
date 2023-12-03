#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/process", showWarnings = FALSE)

parse_size_fnames <- function(x){
 # x <- "data/sponges/data.5000.size_8.gemelli_3.euclidean.permanova.tsv"

  read_tsv(x,
          col_types = cols(.default = col_character())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data_dir", "dataset", "data", "depth", "size_rep",
                      "process", "distance", "dist", "gz")) %>%
    separate(size_rep, into = c("nsamples", "rep")) %>%
    select(dataset, depth, nsamples, rep, process, distance, effect,
          f_statistic, p_value, r2)

}

parse_all_fnames <- function(x){
 # x <- "data/sponges/data.1000.all.raw.bray.permanova.tsv"

  read_tsv(x,
          col_types = cols(.default = col_character())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
             into = c("data_dir", "dataset", "data", "depth", "nsamples",
                      "process", "distance", "dist", "gz")) %>%
    mutate(effect = "treatment)",
           rep = NA_character_) %>%
    select(dataset, depth, nsamples, rep, process, distance, effect,
          f_statistic, p_value, r2)

}

parse_nosize_fnames <- function(x){
 # x <- "data/sponges/data.1000.78_9.sub.bray.permanova.tsv"

  read_tsv(x,
          col_types = cols(.default = col_character())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data_dir", "dataset", "data", "depth",
                      "size_rep", "process", "distance", "dist", "gz")) %>%
    separate(size_rep, into = c("nsamples", "rep")) %>%
    mutate(effect = "treatment") %>%
    select(dataset, depth, nsamples, rep, process, distance, effect,
          f_statistic, p_value, r2)

}



permanova_fnames <- list.files(path = c("data/crc", "data/sponges",
                                        "data/sleep_apnea"),
                               pattern = ".permanova.*.tsv",
                               full.names = TRUE)

size_data <- permanova_fnames %>%
  str_subset(., "size") %>%
  map_dfr(., parse_size_fnames)

all_data <- permanova_fnames %>%
  str_subset(., "all") %>%
  map_dfr(., parse_all_fnames)

nosize_data <- permanova_fnames %>%
  str_subset(., "size|all", negate = TRUE) %>%
  map_dfr(., parse_nosize_fnames)

bind_rows(size_data, all_data, nosize_data) %>%
  write_tsv("data/process/case_studies.permanova.tsv")
