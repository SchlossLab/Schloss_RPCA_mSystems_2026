#!/usr/bin/env Rscript

library(DESeq2)
library(tidyverse)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/positive/data.4000.100_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "deseq2.tsv")

fread(count_fname, colClasses = c(sample = "character")) %>%
  pivot_longer(-sample) %>%
  # filter(value > 0) %>%
  mutate(value = value + 1) %>%  # can't have zero counts
  pivot_wider(names_from = "name", values_from = "value", values_fill = 1) %>%
  column_to_rownames("sample") %>%
  t() %>%
  varianceStabilizingTransformation(fitType = "mean") %>%
  t() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(output_fname)
