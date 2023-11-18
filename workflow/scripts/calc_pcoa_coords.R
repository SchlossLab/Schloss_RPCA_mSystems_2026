#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_file_name <- args[1]
ordination_file_name <- args[2]

# dmatrix_file_name <- "data/simulation/positive_raw.4000.braycurtis.dist.gz"
# ordination_file_name <- "data/simulation/positive_raw.4000.braycurtis.pcoa.axes"

pcoa <- read_tsv(dmatrix_file_name,
                 col_types = cols(sample = col_character(),
                                  .default = col_double())) %>%
  column_to_rownames("sample") %>%
  as.dist() %>%
  cmdscale(eig = TRUE)

points <- pcoa$points

explained <- format(round(100 * pcoa$eig[seq_len(ncol(points))] / sum(pcoa$eig),
                          digits = 2), nsmall = 2L)

colnames(points) <- c(glue("axis1_{explained[1]}"),
                      glue("axis2_{explained[2]}"))

as_tibble(points, rownames = "sample") %>%
  write_tsv(ordination_file_name)
