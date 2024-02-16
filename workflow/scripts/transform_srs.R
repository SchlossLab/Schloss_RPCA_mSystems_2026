#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(SRS)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/positive/data.4000.100_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "srs.tsv")

shared <- fread(count_fname, colClasses = c(sample = "character"))

smallest_group <- shared %>%
  pivot_longer(-sample) %>%
  group_by(sample) %>%
  summarize(n = sum(value)) %>%
  summarize(min = min(n)) %>%
  pull(min)

otu_names <- colnames(shared)[-1]

transformed <- shared %>%
  column_to_rownames("sample") %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin = smallest_group) %>%
  t()

colnames(transformed) <- otu_names

as_tibble(transformed, rownames = "sample") %>%
  write_tsv(output_fname)
