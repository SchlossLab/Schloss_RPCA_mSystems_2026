#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

median_iqr <- function(x) {

    tibble(ymin = quantile(x, prob = 0.25),
            y = median(x),
            ymax = quantile(x, prob = 0.75))
}

datasets <- c("negative", "overlap", "positive")
processes <- c("raw", "oclr", "nclr", "zclr", "gemelli_2", "optspace")
depths <- c(1000, 2000, 4000, 10000)

read_tsv("data/pooled_data/permanova_treatment.tsv") %>%
  select(dataset, process, depth, process, distance, f_statistic, p_value) %>%
  filter(dataset %in% datasets,
        process %in% processes) %>%
  mutate(process_dist = glue("{process}_{distance}"),
        depth = factor(as.character(depth), levels = depths)) %>%
  ggplot(aes(x = depth, y = f_statistic, color = process_dist)) +
  # geom_boxplot() +
  stat_summary(fun.data = median_iqr, size = 0.5) +
  facet_grid(process_dist~dataset, scales = "free")

ggsave("figures/permanova_treatment_sim.tiff",
      width = 7.5, height = 10,
      compression = "lzw+p")
