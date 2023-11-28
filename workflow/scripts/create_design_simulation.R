#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

n_group1 <- as.numeric(args[1])
n_group2 <- as.numeric(args[2])
output_file <- args[3]

n_samples <- n_group1 + n_group2


tibble(index = (1:n_samples) - 1,
       sample = glue("sample_{index}"),
       treatment = c(rep("A", n_group1), rep("B", n_group2))
       ) %>%
  select(-index) %>%
  write_tsv(output_file)
