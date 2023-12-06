#!/usr/bin/env Rscript

library(tidyverse)

pooled_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(dataset %in% c("positive", "overlap", "negative")) %>%
  filter(depth == 4000) %>%
  select(process, distance, dataset, rep, knn_approach, iteration, k_size,
        test, validate) %>%
  mutate(process_calculator = paste(process, distance, sep = "_")) %>%
  filter(process_calculator %in% c("gemelli_2_euclidean", "oclr_euclidean",
                                  "optspace_euclidean", "raw_bray",
                                  "raw_euclidean")) %>%
  filter(rep == 1) %>% #needs to change...
  mutate(process_calculator = factor(process_calculator,
                                     levels = c("optspace_euclidean",
                                                "raw_euclidean",
                                                "oclr_euclidean",
                                                "gemelli_2_euclidean",
                                                "raw_bray")))

format_simulation <- c(
  positive = "Positive",
  overlap = "Overlap",
  negative = "Negative"
  )

pooled_data %>%
  ggplot(aes(x = process_calculator, y = test, fill = knn_approach)) +
  geom_hline(yintercept = 0.5, color = "gray") +
  geom_boxplot() +
  facet_wrap(~dataset, ncol = 1,
             labeller = labeller(dataset = format_simulation)) +
  labs(x = NULL,
       y = "Accuracy") +
  scale_x_discrete(breaks = c("oclr_euclidean",
                              "raw_euclidean", "optspace_euclidean",
                              "raw_bray", "gemelli_2_euclidean"),
                   labels = c("Aitchison", "Raw\nEuclidean",
                              "Raw w/ OptSpace\nEuclidean",
                              "Bray-Curtis", "Robust\nAitchison")) +
  scale_fill_manual(name = NULL,
                    breaks = c("knn_cmfixed", "knn_cmgrid", "knn_pdsgrid"),
                    labels = c("Martino", "Martino+Grid", "Current"),
                    values = c("tomato", "forestgreen", "dodgerblue")) +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"))

ggsave("figures/knn_control_comparison.png", width = 6, height = 5)
