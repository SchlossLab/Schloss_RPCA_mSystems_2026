#!/usr/bin/env Rscript

library(metagenomeSeq)
library(tidyverse)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/positive/data.4000.100_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "css.tsv")

fread(count_fname, colClasses = c(sample = "character")) %>%
  pivot_longer(-sample) %>%
  filter(value > 0) %>%
  pivot_wider(names_from = "name", values_from = "value", values_fill = 0) %>%
  column_to_rownames("sample") %>%
  t() %>%
  newMRexperiment(counts = .) %>% # samples need to be in columns and otus need
  MRcounts(norm = TRUE) %>%       # to be in rows
  t() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(output_fname)
