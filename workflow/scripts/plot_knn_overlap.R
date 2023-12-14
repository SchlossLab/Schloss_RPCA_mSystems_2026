#!/usr/bin/env Rscript

library(tidyverse)

pooled_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(dataset == "overlap") %>%
  mutate(depth = factor(depth, levels = c(1000, 2000, 4000, 10000)),
        process_calculator = paste(process, distance, sep = "_")) %>%
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

format_depth <- c(
  "1000" = "1,000 sequences",
  "2000" = "2,000 sequences",
  "4000" = "4,000 sequences",
  "10000" = "10,000 sequences")

pooled_data %>%
  ggplot(aes(x = process_calculator, y = test, fill = knn_approach)) +
  geom_hline(yintercept = 0.5, color = "gray") +
  geom_boxplot() +
  facet_wrap(~depth, ncol = 1,
             labeller = labeller(depth = format_depth)) +
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
                    labels = c("60/40 Fixed", "60/40 Grid", "80/20 Grid"),
                    values = c("tomato", "forestgreen", "dodgerblue")) +
  scale_y_continuous(limits = c(0.45, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"))

ggsave("figures/knn_overlap_comparison.png", width = 6, height = 8)
