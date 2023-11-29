#!/usr/bin/env Rscript

library(tidyverse)

pooled_data <- read_tsv("data/process/simulation.knn.tsv") %>%
  filter(simulation == "overlap") %>%
  mutate(depth = factor(depth, levels = c(500, 1000, 2000, 4000, 10000)),
         process_calculator = factor(process_calculator,
                                     levels = c("optspace_euclidean",
                                                "raw_euclidean",
                                                "oclr_euclidean",
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
  scale_x_discrete(breaks = c("oclr_euclidean",
                              "raw_euclidean", "optspace_euclidean",
                              "raw_bray", "rpca_euclidean"),
                   labels = c("Aitchison", "Raw\nEuclidean",
                              "Raw w/ OptSpace\nEuclidean",
                              "Bray-Curtis", "Robust\nAitchison")) +
  scale_fill_manual(name = NULL,
                    breaks = c("cm", "cmgrid", "pds"),
                    labels = c("Martino", "Martino+Grid", "Current"),
                    values = c("tomato", "forestgreen", "dodgerblue")) +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"))

ggsave("figures/knn_overlap_comparison.png", width = 6, height = 8)
