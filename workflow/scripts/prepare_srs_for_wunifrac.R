#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ape)
source("workflow/scripts/ape_patch.R")

input <- commandArgs(trailingOnly = TRUE)
count_fname <- input[1]
tree_fname <- input[2]

# count_fname <- "data/crc/data.10439.all_1.srs.tsv.gz"
# tree_fname <- "data/crc/data.10439.all_1.raw.tre"

output_count_fname <- str_replace(count_fname, "srs.tsv.gz", "srs.count_table")


shared <- fread(count_fname, colClasses = c(sample = "character"))


srs_data <- shared %>%
  pivot_longer(-sample, names_to = "Representative_Sequence") %>%
  mutate(total = sum(value), .by = Representative_Sequence) %>%
  filter(total > 0) %>% 
  pivot_wider(values_from = value, names_from = sample, values_fill = 0)


srs_data %>%
  write_tsv(output_count_fname)

output_tree_fname <- str_replace(tree_fname, "raw.tre", "srs.tre")

input_otus <- colnames(shared)[-1]
output_otus <- srs_data %>% pull(Representative_Sequence)
dropped_otus <- input_otus[!input_otus %in% output_otus]

read.tree(tree_fname) %>%
  drop.tip(dropped_otus) %>%
  write.tree(output_tree_fname)
