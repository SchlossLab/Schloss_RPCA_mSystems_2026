#!/usr/bin/env Rscript

library(tidyverse)

pooled_data <- list.files(path = "data/simulation", pattern = "knn_.*.tsv", full.names = TRUE) %>%
  read_tsv(col_types = cols(.default = col_double()),
           id = "file_name") %>%
  mutate(depth = str_replace(file_name,
                             ".*\\.(\\d*)\\..*",
                             "\\1"),
         calculator = str_replace(file_name,
                             ".*\\.\\d*\\.(.*)\\.knn.*",
                             "\\1"),
         process = str_replace(file_name,
                             "data/simulation/overlap_(.*)\\.\\d*\\..*",
                             "\\1"),
         approach = str_replace(file_name,
                                ".*knn_(.*)\\..*",
                              "\\1")
  ) %>%
  select(process, approach, depth, calculator,
         iteration, k_size, test, validate) %>%
  mutate(process_calculator = paste(process, calculator, sep = "_")) %>%
  mutate(depth = factor(depth, levels = c(500, 1000, 2000, 4000, 10000)),
         process_calculator = factor(process_calculator,
                                     levels = c("optspace_euclidean",
                                                "raw_euclidean",
                                                "rpca_euclidean",
                                                "raw_bray")))

format_depth <- c(
  "500" = "500 sequences", 
  "1000" = "1,000 sequences",
  "2000" = "2,000 sequences", 
  "4000" = "4,000 sequences",
  "10000" = "10,000 sequences")

pooled_data %>%
  ggplot(aes(x = process_calculator, y = test, fill = approach)) +
  geom_hline(yintercept = 0.5, color = "gray") +
  geom_boxplot() +
  facet_wrap(~depth, ncol = 1,
             labeller = labeller(depth = format_depth)) +
  labs(x = NULL,
       y = "Accuracy") +
  scale_x_discrete(breaks = c("raw_euclidean", "optspace_euclidean",
                              "raw_bray", "rpca_euclidean"),
                   labels = c("Raw\nEuclidean", "Raw w/ OptSpace\nEuclidean",
                              "Raw\nBray", "RCLR w/ OptSpace\nEuclidean")) +
  scale_fill_manual(name = NULL,
                    breaks = c("cm", "pds"),
                    labels = c("Martino", "Current"),
                    values = c("hotpink", "dodgerblue")) +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"))

ggsave("figures/overlap_knn_analysis.png", width = 5, height = 8)
