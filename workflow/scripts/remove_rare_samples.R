#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ape)
library(glue)
source("workflow/scripts/ape_patch.R")

args <-  commandArgs(trailingOnly = TRUE)
count_meta_fname <- args[1]
tree_fname <- args[2]
n_sequence_threshold <- as.numeric(args[3])

# count_meta_fname <- "data/crc/data_treatment.tsv.gz"
# tree_fname <- "data/crc/data_treatment.tre"
# n_sequence_threshold <- 10000

count_fname_out <- str_replace(
  count_meta_fname,
  "_treatment.*",
  glue(".{n_sequence_threshold}.all_1.raw.tsv.gz"))

tree_fname_out <- str_replace(
  count_meta_fname,
  "_treatment.*",
  glue(".{n_sequence_threshold}.all_1.raw.tre"))

counttable_fname_out <- str_replace(
  count_meta_fname,
  "_treatment.*",
  glue(".{n_sequence_threshold}.all_1.raw.count_table"))

design_fname_out <- str_replace(
  count_meta_fname,
  "_treatment.*",
  glue(".{n_sequence_threshold}.all_1.design"))


count_meta <- fread(count_meta_fname)

filtered_data <- count_meta %>%
  pivot_longer(cols = -c(sample, treatment)) %>%
  mutate(n_seqs = sum(value), .by = "sample") %>%
  filter(n_seqs >= n_sequence_threshold) %>%
  filter(value > 0) %>%
  select(-n_seqs) %>%
  pivot_wider(id_cols = c(sample, treatment),
              names_from = name,
              values_from = value, values_fill = 0)

otu_names_in <- colnames(count_meta)[-c(1, 2)]
otu_names_out <- colnames(filtered_data)[-c(1, 2)]
otu_names_remove <- otu_names_in[!otu_names_in %in% otu_names_out]

filtered_data %>%
  select(-treatment) %>%
  write_tsv(count_fname_out)

filtered_data %>%
  select(sample, treatment) %>%
  write_tsv(design_fname_out)

tr <- read.tree(tree_fname)

# for some reason clearcut has a node without a branch length. i assume this is
# because it is unrooted and should be assumed to have a zero branch length
tr$edge.length[is.nan(tr$edge.length)] <- 0

tr %>%
  drop.tip(otu_names_remove) %>%
  write.tree_patch(tree_fname_out)

filtered_data %>%
  select(-treatment) %>%
  pivot_longer(-sample, names_to = "Representative_Sequence") %>%
  mutate(total = sum(value), .by = Representative_Sequence) %>%
  pivot_wider(values_from = value, names_from = sample) %>%
  write_tsv(counttable_fname_out)
