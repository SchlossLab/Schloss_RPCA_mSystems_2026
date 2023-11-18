#!/usr/bin/env Rscript

library(tidyverse)
library(broom)

args <-  commandArgs(trailingOnly = TRUE)

csv_file_name <- args[1]
output_file_name <- args[2]

read_csv(csv_file_name) %>%
  select(sequence_depth, OTUs, starts_with("sample")) %>%
  pivot_longer(cols = starts_with("sample")) %>%
  group_by(sequence_depth, name) %>%
  summarize(n = sum(value), .groups = "drop") %>%
  mutate(name = str_replace(name, "sample_", "") %>% as.numeric) %>%
  arrange(sequence_depth, name) %>%
  mutate(treatment = rep(rep(c("A", "B"), each = 100), 5)) %>%
  nest(data = -sequence_depth) %>%
  mutate(test = map(data, ~wilcox.test(.x$n ~ .x$treatment) %>% tidy)) %>%
  select(sequence_depth, test) %>%
  unnest(test) %>%
  select(sequence_depth, p_value = p.value) %>%
  write_tsv(output_file_name)
