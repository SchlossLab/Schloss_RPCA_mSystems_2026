#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

median_iqr <- function(x) {

    tibble(ymin = quantile(x, prob = 0.25),
            y = median(x),
            ymax = quantile(x, prob = 0.75))
}

datasets <- c("sponges", "sleep_apnea", "crc")
processes <- c("raw", "sub", "rare", "oclr", "nclr", "zclr", "optspace")

results_summary <- read_tsv("data/pooled_data/permanova_treatment.tsv") %>%
  filter(nsamples_per_group == "all",
        dataset %in% datasets,
        process %in% processes | str_detect(process, "gemelli_")) %>%
  mutate(depth_chr = if_else(depth == 1000, "low", "high")) %>%
  select(dataset, process, distance, depth_chr, p_value, f_statistic)



results_summary %>%
  filter(dataset != 'crc') %>%
  mutate(test = p_value < 0.05) %>%
  pull(test) %>%
  stopifnot(all(.))


results_summary %>%
  filter(dataset == 'crc') %>%
  select(-f_statistic) %>%
  pivot_wider(values_from = p_value, names_from = depth_chr) %>%
  filter(low > 0.05 | high > 0.05)

results_summary %>%
  select(-p_value) %>%
  pivot_wider(names_from = "dataset", values_from = "f_statistic") %>%
  print(n = Inf)
