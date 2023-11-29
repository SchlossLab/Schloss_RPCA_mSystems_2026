#!/usr/bin/env Rscript

library(tidyverse)

pooled_data <- read_tsv("data/process/simulation.knn.tsv") %>%
#  filter(simulation != "overlap") %>%
  filter(depth == 4000) %>%
  select(process_calculator, simulation, approach, iteration, k_size, test, validate) %>%
  mutate(process_calculator = factor(process_calculator,
                                     levels = c("optspace_euclidean",
                                                "raw_euclidean",
                                                "oclr_euclidean",
                                                "rpca_euclidean",
                                                "raw_bray")))

format_simulation <- c(
  positive = "Positive",
  overlap = "Overlap",
  negative = "Negative"
  )

pooled_data %>%
  ggplot(aes(x = process_calculator, y = test, fill = approach)) +
  geom_hline(yintercept = 0.5, color = "gray") +
  geom_boxplot() +
  facet_wrap(~simulation, ncol = 1,
             labeller = labeller(simulation = format_simulation)) +
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

ggsave("figures/knn_control_comparison.png", width = 6, height = 5)
