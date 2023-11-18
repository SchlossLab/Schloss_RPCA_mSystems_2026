#!/usr/bin/env Rscript

library(tidyverse)

list.files(path = "data/simulation", pattern = "knn.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(distance_file = col_character(),
                            replicate = col_integer(),
                            accuracy = col_double())) %>%
  mutate(depth = str_replace(distance_file,
                             ".*\\.(\\d*)\\..*",
                             "\\1"),
         calculator = str_replace(distance_file,
                             ".*\\.\\d*\\.(.*).dist.gz",
                             "\\1"),
         process = str_replace(distance_file,
                             "data/simulation/positive_(.*)\\.(\\d*)\\..*",
                             "\\1")
  ) %>%
  select(process, depth, calculator, replicate, accuracy) %>% 
  mutate(process_calculator = paste(process, calculator, sep = "_")) %>%
  mutate(depth = factor(depth, levels = c(500, 1000, 2000, 4000, 10000))) %>%
  ggplot(aes(x = process_calculator, y = accuracy,
             fill = depth, color = depth)) +
  geom_boxplot()

ggsave("figures/positive_knn_analysis.png", width = 7, height = 4)