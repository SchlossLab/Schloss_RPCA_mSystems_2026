#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

pcoa_fname <- args[1]

# pcoa_fname <- "data/overlap/data.4000.100_1.rare.bray.pcoa.axes"

output_fname <- str_replace(pcoa_fname,
                            "rare.bray.pcoa.axes",
                            "bray_pcoa.euclidean.dist.gz")

read_tsv(pcoa_fname) %>%
  select(sample, starts_with("axis1_"), starts_with("axis2_")) %>%
  column_to_rownames("sample") %>%
  vegdist(method = "euclidean") %>%
  as.matrix() %>%
  as.data.frame() %>%
  as_tibble(rownames = "sample") %>%
  write_tsv(output_fname)
