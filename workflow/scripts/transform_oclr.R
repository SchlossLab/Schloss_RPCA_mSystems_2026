#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/positive/data.4000.100_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "oclr.tsv")



geometric_mean <- function(x) {

  exp(mean(log(x[x > 0])))

}

fread(count_fname, colClasses = c(sample = "character")) %>%
  pivot_longer(-sample) %>%
    group_by(sample) %>%
    mutate(value = value + 1,
           oclr = log(value / geometric_mean(value))) %>%
    ungroup() %>%
    select(-value) %>%
    pivot_wider(names_from = name, values_from = oclr, values_fill = 0) %>%
    write_tsv(output_fname)
