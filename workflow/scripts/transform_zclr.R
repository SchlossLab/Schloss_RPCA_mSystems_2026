#!/usr/bin/env Rscript

library(zCompositions) #has to go first since it has a different select than tdv
library(tidyverse)
library(data.table)

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]

# count_fname <- "data/crc/data.10439.all_1.raw.tsv.gz"

output_fname <- str_replace(count_fname, "raw.tsv", "zclr.tsv")



geometric_mean <- function(x) {

  exp(mean(log(x[x > 0])))

}

fread(count_fname, colClasses = c(sample = "character")) %>%
  column_to_rownames("sample") %>%
  cmultRepl(method = "CZM", output = "p-count") %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  mutate(zclr = log(value / geometric_mean(value)), .by = sample) %>%
  select(-value) %>%
  pivot_wider(names_from = name, values_from = zclr, values_fill = 0) %>%
  write_tsv(output_fname)
