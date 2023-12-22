#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(ggh4x)

args <- commandArgs(trailingOnly = TRUE)
fitting_method <- args[1]

median_iqr <- function(x) {

    tibble(ymin = quantile(x, prob = 0.25),
            y = median(x),
            ymax = quantile(x, prob = 0.75))
}


processes_distances <- c(
  "raw.euclidean" = "Euclidean\n(Raw)",
  "optspace.euclidean" = "OptSpace",
  "oclr.euclidean" = "One CLR",
  "nclr.euclidean" = "Nudge CLR",
  "zclr.euclidean" = "Zero CLR",
  "gemelli_2.euclidean" = "RPCA\n(2 parts)",
  "gemelli_auto.euclidean" = "RPCA\n(auto parts)",
  "rare.bray" = "Bray-Curtis\n(Rarefaction)",
  "raw.bray" = "Bray-Curtis\n(Raw)",
  "sub.bray" = "Bray-Curtis\n(Subsampled)"
)

case_studies <- c(
  sponges = "Sponges",
  sleep_apnea = "Sleep apnea",
  crc = "CRC"
)

depths <- c(
  low = "Low",
  high = "High"
)



pooled_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  mutate(process_distance = glue("{process}.{distance}"),
        depth_chr = if_else(depth == 1000, "low", "high")) %>%
  filter(knn_approach == fitting_method,
         dataset %in% names(case_studies),
         process_distance %in% names(processes_distances)) %>%
  select(dataset, depth_chr, nsamples_per_group, process_distance,
        rep, iteration, k_size, test, validate) %>%
  mutate(process_distance = factor(process_distance,
                                   levels = names(processes_distances),
                                   labels = processes_distances),
        depth_chr = factor(depth_chr,
                           levels = names(depths),
                           labels = depths),
        dataset = factor(dataset,
                         levels = names(case_studies),
                         labels = case_studies))#,
        #nsamples_per_group = factor(nsamples_per_group))

pooled_data %>%
  ggplot(aes(x = nsamples_per_group, y = test,
            group = glue("{rep}_{depth_chr}"),
            color = depth_chr)) +
  geom_hline(yintercept = 0.5, color = "gray") +
  stat_summary(fun.data = median_iqr, size = 0.1,
              position = position_dodge(width = 0.4), alpha = 0.25) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
              group = 1, geom = "crossbar", width = 1,
              show.legend = FALSE) +
  facet_grid2(rows = vars(dataset), cols = vars(process_distance),
              scales = "free_x", independent = "x") +
  labs(x = "Number of samples per group",
       y = "Accuracy") +
  theme_classic() +
  theme(strip.background = element_rect(color = NA, fill = "lightgray"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.1),
        legend.position = "bottom")

ggsave(glue("figures/{fitting_method}_cs.tiff"),
      width = 9, height = 5,
      compression = "lzw+p")
