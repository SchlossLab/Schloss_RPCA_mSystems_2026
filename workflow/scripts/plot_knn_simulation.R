#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
fitting_method <- args[1]

median_iqr <- function(x) {

    tibble(ymin = quantile(x, prob = 0.25),
            y = median(x),
            ymax = quantile(x, prob = 0.75))
}


processes_distances <- c(
  "raw.euclidean" = "Euclidean",
  "optspace.euclidean" = "OptSpace",
  "oclr.euclidean" = "One CLR",
  "nclr.euclidean" = "Nudge CLR",
  "zclr.euclidean" = "Zero CLR",
  "gemelli_2.euclidean" = "RPCA",
  "raw.bray" = "Bray-Curtis"
)

simulations <- c(
  positive = "Positive",
  overlap = "Overlap",
  negative = "Negative"
)

depths <- c(
  "1000" = "1,000",
  "2000" = "2,000",
  "4000" = "4,000",
  "10000" = "10,000"
)



pooled_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  mutate(process_distance = glue("{process}.{distance}")) %>%
  filter(knn_approach == fitting_method,
         dataset %in% names(simulations),
         process_distance %in% names(processes_distances)) %>%
  select(dataset, depth, process_distance, rep, iteration,
        k_size, test, validate) %>%
  mutate(process_distance = factor(process_distance,
                                   levels = names(processes_distances),
                                   labels = processes_distances),
        depth = factor(as.character(depth),
                        levels = names(depths),
                        labels = depths),
        dataset = factor(dataset,
                         levels = names(simulations),
                         labels = simulations))

pooled_data %>%
  ggplot(aes(x = depth, y = test, group = glue("{rep}_{depth}"))) +
  geom_hline(yintercept = 0.5, color = "gray") +
  stat_summary(fun.data = median_iqr, size = 0.25, color = "gray",
              position = position_dodge(width = 0.4)) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
                color = "red", group = 1, geom = "crossbar", width = 1) +
  facet_grid(dataset ~ process_distance) +
  labs(x = "Number of sequences",
       y = "Accuracy") +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.1))

ggsave(glue("figures/{fitting_method}_sim.tiff"),
      width = 7, height = 5,
      compression = "lzw+p")
