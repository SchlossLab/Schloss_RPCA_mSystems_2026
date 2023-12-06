#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)

args <-  commandArgs(trailingOnly = TRUE)
count_in_fname <- args[1]
design_in_fname <- args[2]
nsamples <- as.integer(args[3])
rep <- as.integer(args[4])

# count_in_fname <- "data/crc/data.1000.all_1.raw.tsv.gz"
# design_in_fname <- "data/crc/data.1000.all_1.design"
# nsamples <- 35
# rep <- 1

set.seed(rep)

count_out_fname <- str_replace(count_in_fname,
                               "all_1.raw.tsv.gz",
                               glue("{nsamples}_{rep}.raw.tsv.gz"))

design_out_fname <- str_replace(count_in_fname,
                                "all_1.raw.tsv.gz",
                                glue("{nsamples}_{rep}.design"))

count_data <- fread(count_in_fname,
                    colClasses = c("sample" = "character"))

design_data <- read_tsv(design_in_fname,
                        col_types = cols(sample = col_character(),
                                        treatment = col_character()))

count_design <- inner_join(design_data, count_data, by = "sample") %>%
  slice_sample(n = nsamples, by = treatment) %>%
  pivot_longer(cols = starts_with("OTU")) %>%
  filter(value > 0) %>%
  pivot_wider(id_cols  = c(sample, treatment), values_fill = 0)

count_design %>%
  select(sample, starts_with("OTU")) %>%
  write_tsv(count_out_fname)

count_design %>%
  select(sample, treatment) %>%
  write_tsv(design_out_fname)
