#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/positive/data.4000.100_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "relabund.tsv")

fread(count_fname, colClasses = c(sample = "character")) %>%
  pivot_longer(-sample) %>%
  mutate(relabund = value / sum(value), .by = sample) %>%
  select(-value) %>%
  pivot_wider(values_from = relabund, names_from = name, values_fill = 0) %>%
  write_tsv(output_fname)
