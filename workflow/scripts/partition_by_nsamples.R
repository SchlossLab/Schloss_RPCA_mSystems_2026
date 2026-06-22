#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(glue)
library(ape)
source("workflow/scripts/ape_patch.R")

args <-  commandArgs(trailingOnly = TRUE)
count_in_fname <- args[1]
tree_in_fname <- args[2]
design_in_fname <- args[3]
nsamples <- as.integer(args[4])
rep <- as.integer(args[5])

# count_in_fname <- "data/sleep_apnea/data.11879.all_1.raw.tsv.gz"
# tree_in_fname <- "data/sleep_apnea/data.11879.all_1.raw.tre"
# design_in_fname <- "data/sleep_apnea/data.11879.all_1.design"
# nsamples <- 35
# rep <- 10

set.seed(rep)

count_out_fname <- str_replace(count_in_fname,
                               "all_1.raw.tsv.gz",
                               glue("{nsamples}_{rep}.raw.tsv.gz"))

tree_out_fname <- str_replace(tree_in_fname,
                              "all_1.raw.tre",
                              glue("{nsamples}_{rep}.raw.tre"))

counttable_out_fname <- str_replace(tree_in_fname,
                                    "all_1.raw.tre",
                                    glue("{nsamples}_{rep}.raw.count_table"))

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
  pivot_longer(cols = -c(sample, treatment)) %>%
  filter(value > 0) %>%
  pivot_wider(id_cols  = c(sample, treatment), values_fill = 0)

count_design %>%
  select(-treatment) %>%
  write_tsv(count_out_fname)

count_design %>%
  select(sample, treatment) %>%
  write_tsv(design_out_fname)

in_otu_names <- colnames(count_data)[-1]
out_otu_names <- colnames(count_design)[-c(1, 2)]
remove_otu_names <- in_otu_names[!in_otu_names %in% out_otu_names]

options(scipen = 9999999)

read.tree(tree_in_fname) %>%
  drop.tip(remove_otu_names) %>%
  write.tree_patch(tree_out_fname)

count_design %>%
  select(-treatment) %>%
  pivot_longer(-sample, names_to = "Representative_Sequence") %>%
  mutate(total = sum(value), .by = Representative_Sequence) %>%
  pivot_wider(values_from = value, names_from = sample) %>%
  write_tsv(counttable_out_fname)
