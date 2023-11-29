#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

args <-  commandArgs(trailingOnly = TRUE)
count_meta_fname <- args[1]
nsamples_rep <- args[2]

# count_meta_fname <- "data/sponges/data.13755.raw_design.tsv.gz"
# nsamples_rep <- "55_2"

nsamples_rep_parse <- as.numeric(unlist(str_split(nsamples_rep, pattern = "_")))
n_samples <- nsamples_rep_parse[1]
rep <- nsamples_rep_parse[2]

count_fname <- str_replace(count_meta_fname,
                          "raw_design.tsv.gz",
                          glue("{nsamples_rep}.raw.tsv.gz"))

design_fname <- str_replace(count_meta_fname,
                            "raw_design.tsv.gz",
                            glue("{nsamples_rep}.design"))


set.seed(rep)

count_design <- fread(count_meta_fname) %>%
  slice_sample(n = n_samples, by = treatment) %>%
  pivot_longer(cols = starts_with("OTU")) %>%
  filter(value > 0) %>%
  pivot_wider(id_cols  = c(sample, treatment), values_fill = 0)

count_design %>%
  select(sample, starts_with("OTU")) %>%
  write_tsv(count_fname)

count_design %>%
  select(sample, treatment) %>%
  write_tsv(design_fname)
