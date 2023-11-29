#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(glue)

set.seed(19760620)

args <- commandArgs(trailingOnly = TRUE)

count_table_name <- args[1]
distance_method <- args[2]
seq_depth <- as.numeric(args[3])

# count_table_name <- "data/sponges/data.1000.78_9.raw.tsv.gz"
# distance_method <- "bray"
# seq_depth <- 1000

output_fname <- str_replace(count_table_name,
                            "tsv.gz",
                            glue("sub_{distance_method}.dist.gz"))

read_tsv(count_table_name) %>%
  column_to_rownames("sample") %>%
  avgdist(sample = seq_depth, dmethod = distance_method, iterations = 1) %>%
  as.matrix() %>%
  as.data.frame() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(output_fname)
