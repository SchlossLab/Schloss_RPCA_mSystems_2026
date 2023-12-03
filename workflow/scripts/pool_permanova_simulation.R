#!/usr/bin/env Rscript

library(tidyverse)

dir.create("data/process", showWarnings = FALSE)


parse_size_fnames <- function(x){
 # x <- "data/simulation/positive.5000.size_8.optspace.euclidean.permanova.tsv"

  read_tsv(x,
          col_types = cols(dmatrix_fname = col_character(),
                          effect = col_character(),
                          .default = col_double())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data", "simulation", "dataset", "depth", "size_rep",
                      "process", "distance", "dist", "gz")) %>%
    select(dataset, depth, size_rep, process, distance, effect, f_statistic,
          p_value, r2) %>%
    mutate(size_rep = str_replace(size_rep, "size_", "") %>% as.numeric)

}

parse_nosize_fnames <- function(x){
 # x <- "data/simulation/overlap.4000.optspace.euclidean.permanova.tsv"

  read_tsv(x,
          col_types = cols(dmatrix_fname = col_character(),
                          .default = col_double())) %>%
    separate(dmatrix_fname, sep = "[\\/\\.]",
              into = c("data", "simulation", "dataset", "depth",
                      "process", "distance", "dist", "gz")) %>%
    mutate(size_rep = NA_real_,
          effect = "treatment") %>%
    select(dataset, depth, size_rep, process, distance, effect, f_statistic,
          p_value, r2)

}



permanova_fnames = list.files(path = "data/simulation", pattern = ".permanova.*.tsv", full.names = TRUE)

size_data <- permanova_fnames %>%
  str_subset(., "size") %>%
  map_dfr(., parse_size_fnames)

nosize_data <- permanova_fnames %>%
  str_subset(., "size", negate = TRUE) %>%
  map_dfr(., parse_nosize_fnames)

bind_rows(size_data, nosize_data) %>%
  write_tsv("data/process/simulation.permanova.tsv")


# size_data %>%
#   select(-r2, -f_statistic) %>%
#   pivot_wider(names_from = effect, values_from = "p_value") %>%
#   group_by(dataset, depth, process, distance) %>%
#   summarize(treatment = mean(treatment < 0.05),
#             size = mean(size < 0.05),
#             interaction = mean(interaction < 0.05), .groups = "drop")
