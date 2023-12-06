#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/pooled_data", showWarnings = FALSE)

pcoa_axes_fnames <- commandArgs(trailingOnly = TRUE)

# pcoa_axes_fnames <- c(
#   "data/crc/data.10439.83_1.gemelli_auto.euclidean.pcoa.axes",
#   "data/crc/data.10439.83_2.gemelli_auto.euclidean.pcoa.axes",
#   "data/crc/data.10439.83_3.gemelli_auto.euclidean.pcoa.axes",
#   "data/crc/data.10439.83_4.gemelli_auto.euclidean.pcoa.axes")


process_file <- function(x){

  split_fname <- str_split(x, "[\\.\\/]") %>% unlist()
  n_axes <- read_tsv(x) %>% ncol() - 1

  tibble(
      dataset = split_fname[2],
      depth =  split_fname[4],
      samples_rep = split_fname[5],
      distance = split_fname[7],
      n_axes = n_axes
    ) %>%
    separate(samples_rep, into = c("nsamples_per_group", "rep"))

}

map_dfr(pcoa_axes_fnames, process_file) %>%
  write_tsv("data/pooled_data/gemelli_auto_naxes.tsv")
