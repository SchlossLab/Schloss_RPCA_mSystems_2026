#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_file_name <- args[1]
design_file_name <- args[2]
output_file_name <- args[3]

# dmatrix_file_name <- "data/simulation/overlap_raw.4000.braycurtis.dist.gz"
# design_file_name <- "data/simulation/overlap_raw.design.gz"

dmatrix <- read_tsv(dmatrix_file_name,
                    col_types = cols(sample = col_character(),
                                     .default = col_double())) %>%
  column_to_rownames("sample") %>%
  as.dist()

design <- read_tsv(design_file_name,
                   col_types = cols(.default = col_character())) %>%
  column_to_rownames("sample")

stopifnot(all(rownames(design) == labels(dmatrix)))

test_result <- adonis2(dmatrix ~ design$treatment)

tibble(dmatrix_file_name = dmatrix_file_name,
       f_statistic = test_result$F[1],
       p_value = test_result$`Pr(>F)`[1],
       r2 = test_result$R2[1]) %>%
  write_tsv(output_file_name)
