#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ape)
library(glue)
source("workflow/scripts/ape_patch.R")

args <-  commandArgs(trailingOnly = TRUE)
count_in_fname <- args[1]
tree_in_fname <- args[2] # if == "no_tree" skips tree processing
design_in_fname <- args[3]
seed <- args[4]

# count_in_fname <- "data/sponges/data.10000.all_1.raw.tsv.gz"
# tree_in_fname <- "data/sponges/data.10000.all_1.raw.tre"
# design_in_fname <- "data/sponges/data.10000.all_1.design"
# seed <- 1

# count_in_fname <- "data/positive/data.10000.100_1.raw.tsv.gz"
# tree_in_fname <- "no_tree"
# design_in_fname <- "data/positive/data.10000.100_1.design"
# seed <- 1


set.seed(seed)

count_out_fname <- str_replace(count_in_fname,
                              "\\.10000\\..*\\.raw.tsv.gz",
                              glue(".5000.size_{seed}.raw.tsv.gz"))

design_out_fname <- str_replace(design_in_fname,
                                "\\.10000\\..*\\.design",
                                glue(".5000.size_{seed}.design"))

count_data <- fread(count_in_fname,
                    colClasses = list(character = "sample"))

design_data <- read_tsv(design_in_fname,
                        col_types = cols(sample = col_character()))


count_design <- inner_join(count_data, design_data, by = "sample") %>%
  pivot_longer(cols = -c(sample, treatment)) %>%
  filter(value > 0) %>%
  uncount(weights = value)

count_design_10k <- count_design %>%
  slice_sample(n = 10000, by = c(sample, treatment), replace = FALSE) %>%
  count(sample, treatment, name) %>%
  mutate(sample = paste0(sample, "_10k"))

count_design_5k <- count_design %>%
  slice_sample(n = 5000, by = c(sample, treatment), replace = FALSE) %>%
  count(sample, treatment, name) %>%
  mutate(sample = paste0(sample, "_5k"))

count_design_pool <- bind_rows(count_design_10k, count_design_5k) %>%
  pivot_wider(names_from = name, values_from = n, values_fill = 0)


count_design_pool %>%
  select(-treatment) %>%
  write_tsv(count_out_fname)

count_design_pool %>%
  select(sample, treatment) %>%
  write_tsv(design_out_fname)

if(tree_in_fname != "no_tree"){

  tree_out_fname <- str_replace(tree_in_fname,
                                "\\.10000\\..*\\.raw.tre",
                                glue(".5000.size_{seed}.raw.tre"))

  counttable_out_fname <- str_replace(tree_in_fname,
                                "\\.10000\\..*\\.raw.tre",
                                glue(".5000.size_{seed}.raw.count_table"))

  in_otu_names <- colnames(count_data)[-1]
  out_otu_names <- colnames(count_design_pool)[-c(1, 2)]
  remove_otu_names <- in_otu_names[!in_otu_names %in% out_otu_names]

  read.tree(tree_in_fname) %>%
    drop.tip(remove_otu_names) %>%
    write.tree_patch(tree_out_fname)

  count_design_pool %>%
    select(-treatment) %>%
    pivot_longer(-sample, names_to = "Representative_Sequence") %>%
    mutate(total = sum(value), .by = Representative_Sequence) %>%
    pivot_wider(values_from = value, names_from = sample) %>%
    write_tsv(counttable_out_fname)

}
