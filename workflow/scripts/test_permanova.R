#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_fname <- args[1]
design_fname <- args[2]
output_fname <- args[3]

# dmatrix_fname <- "data/simulation/negative.4000.optspace.euclidean.dist.gz"
# design_fname <- "data/simulation/simulation.design"

dmatrix_design <- inner_join(
                    read_tsv(dmatrix_fname,
                             col_types = cols(sample = col_character(),
                                              .default = col_double())),
                    read_tsv(design_fname,
                            col_types = cols(.default = col_character())),
                    by = "sample") %>%
                  select(sample, treatment, everything())


dmatrix <- dmatrix_design %>%
  select(-treatment) %>%
  column_to_rownames("sample") %>%
  as.dist()

design <- dmatrix_design %>%
  select(sample, treatment) %>%
  column_to_rownames("sample")

stopifnot(all(rownames(design) == labels(dmatrix)))

test_result <- adonis2(dmatrix ~ design$treatment)

tibble(dmatrix_fname = dmatrix_fname,
       f_statistic = test_result$F[1],
       p_value = test_result$`Pr(>F)`[1],
       r2 = test_result$R2[1]) %>%
  write_tsv(output_fname)
