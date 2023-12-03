#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_fname <- args[1]
design_fname <- args[2]
output_fname <- args[3]

# dmatrix_fname <- "data/sponges/data.5000.size_1.gemelli_2.euclidean.dist.gz"
# design_fname <- "data/sponges/data.5000.size_1.design"

dmatrix_design <- inner_join(
                    read_tsv(dmatrix_fname,
                             col_types = cols(sample = col_character(),
                                              .default = col_double())),
                    read_tsv(design_fname,
                            col_types = cols(.default = col_character())),
                    by = "sample") %>%
                  mutate(size = str_replace(sample, ".*_(\\d*k)", "\\1")) %>%
                  select(sample, treatment, size, everything())


dmatrix <- dmatrix_design %>%
  select(-treatment, -size) %>%
  column_to_rownames("sample") %>%
  as.dist()

design <- dmatrix_design %>%
  select(sample, treatment, size) %>%
  column_to_rownames("sample")

stopifnot(all(rownames(design) == labels(dmatrix)))

test_result <- adonis2(dmatrix ~ design$treatment * design$size)

# NOTE: We don't really care about the treatment effect in this simulation. It's
# p-value will be inflated since we're artificially duplicating every sample.
# This will make each treatment look more similar to itself than it really is,
# which will reduce the intra-treatment variation increasing the F statistic
# and decreasing the p-value for the treatment effect.

tibble(dmatrix_fname = dmatrix_fname,
        effect = c("treatment", "size", "interaction"),
        f_statistic = test_result$F[1:3],
        p_value = test_result$`Pr(>F)`[1:3],
       r2 = test_result$R2[1:3]) %>%
  write_tsv(output_fname)
