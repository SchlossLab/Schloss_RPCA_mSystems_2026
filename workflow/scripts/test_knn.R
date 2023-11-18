#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

dmatrix_file_name <- args[1]
design_file_name <- args[2]
output_file_name <- args[3]

# dmatrix_file_name <- "data/simulation/positive_raw.500.euclidean.dist.gz"
# design_file_name <- "data/simulation/positive_raw.design.gz"
# output_file_name <- "data/simulation/positive_raw.500.euclidean.tsv"

set.seed(19760620)

f_training <- 0.6
f_testing <- 0.4
n_neighbors <- 13

dmatrix <- read_tsv(dmatrix_file_name)
design <- read_tsv(design_file_name)
full_design_dmatrix <- inner_join(design, dmatrix, by = "sample")


run_knn <- function(x) {

# The next three code chunks partitions the data to have an equal number of A
# and B samples in the training and test dataset
partitions <- x %>%
  select(sample, treatment) %>%
  mutate(train_test = sample(c(rep("train", f_training * n()),
                               rep("test", f_testing * n()))),
         .by = "treatment")

train_samples <- partitions %>%
  filter(train_test == "train") %>%
  pull(sample)

test_samples <- partitions %>%
  filter(!sample %in% train_samples) %>%
  pull(sample)

x %>%
  filter(sample %in% train_samples) %>%
  select(train_sample = sample, treatment, all_of(test_samples)) %>%
  pivot_longer(cols = starts_with("sample_"),
               names_to = "test_sample",
               values_to = "distance") %>%
  nest(data = -test_sample) %>%
  mutate(classification = map_chr(data, ~{
      .x %>%
        slice_min(n = n_neighbors, order_by = distance) %>%
        count(treatment) %>%
        mutate(proportion = n / sum(n)) %>%
        slice_max(order_by = proportion, n = 1) %>%
        pull(treatment) %>%
        sample(size = 1) # if there's a tie, randomly pick one treatment group
    }
  )) %>%
  select(test_sample, classification) %>%
  inner_join(., x %>% select(sample, treatment), by = c("test_sample" = "sample")) %>%
  rename(truth = treatment) %>%
  count(classification, truth) %>%
  mutate(correct = classification == truth) %>%
  group_by(correct) %>%
  summarize(n = sum(n)) %>%
  mutate(fraction = n / sum(n)) %>%
  filter(correct == TRUE) %>%
  select(accuracy = fraction)

}

map_dfr(1:100, ~run_knn(full_design_dmatrix), .id = "replicate") %>%
  mutate(distance_file = dmatrix_file_name) %>%
  select(distance_file, replicate, accuracy) %>%
  write_tsv(output_file_name)
