#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_fname <- args[1]
ordination_fname <- args[2]

# dmatrix_fname <- "data/simulation/overlap.4000.raw.braycurtis.dist.gz"
# ordination_fname <- "data/simulation/overlap.4000.raw.braycurtis.pcoa.axes"

pcoa <- read_tsv(dmatrix_fname,
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
  write_tsv(ordination_fname)
