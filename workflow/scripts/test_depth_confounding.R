#!/usr/bin/env Rscript

library(tidyverse)
library(broom)

args <-  commandArgs(trailingOnly = TRUE)
depth_fname <- args[1]
output_fname <- args[2]

# depth_fname <- "data/sponges/data.depths.tsv.gz"
# output_fname <- "data/sponges/data.depth_confounding.tsv"

get_p_value <- function(depth_treats, threshold){

  depth_treats %>%
    filter(n_seqs > threshold) %>%
    wilcox.test(n_seqs ~ treatment, data = .) %>%
    tidy() %>%
    rename(p_value = p.value)

}


depths <- read_tsv(depth_fname, col_types = cols(sample = col_character(),
                                                 treatment = col_character(),
                                                 n_seqs = col_double())) %>%
  arrange(n_seqs)

map_dfr(depths$n_seqs[1:as.integer(nrow(depths) / 2)],
        ~get_p_value(depths, threshold = .x)) %>%
  mutate(n_seqs = depths$n_seqs[1:as.integer(nrow(depths) / 2)],
         .before = statistic) %>%
  write_tsv(output_fname)
