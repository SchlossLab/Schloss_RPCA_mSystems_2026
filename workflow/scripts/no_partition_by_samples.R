#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

args <-  commandArgs(trailingOnly = TRUE)
count_meta_fname <- args[1]

# count_meta_fname <- "data/sponges/data.13755.raw_design.tsv.gz"
# nsamples_rep <- "all"

count_fname <- str_replace(count_meta_fname,
                          "raw_design.tsv.gz",
                          "all.raw.tsv.gz")

design_fname <- str_replace(count_meta_fname,
                            "raw_design.tsv.gz",
                            "all.design")

count_design <- fread(count_meta_fname) %>%
  pivot_longer(cols = starts_with("OTU")) %>%
  filter(value > 0) %>%
  pivot_wider(id_cols  = c(sample, treatment), values_fill = 0)

count_design %>%
  select(sample, starts_with("OTU")) %>%
  write_tsv(count_fname)

count_design %>% 
  select(sample, treatment) %>%
  write_tsv(design_fname)
