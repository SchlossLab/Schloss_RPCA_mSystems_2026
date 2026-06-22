#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(glue)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

count_table_name <- args[1]
distance_method <- args[2]

# count_table_name <- "data/simulation/overlap_raw.500.tsv.gz"
# distance_method <- "bray"

output_fname <- str_replace(count_table_name,
                                "tsv.gz",
                                glue("{distance_method}.dist.gz"))

fread(count_table_name) %>%
  column_to_rownames("sample") %>%
  vegdist(method = distance_method) %>%
  as.matrix() %>%
  as.data.frame() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(output_fname)
