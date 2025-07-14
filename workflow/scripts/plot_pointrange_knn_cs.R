#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(ggtext)

pretty_processes <- c(
  "raw" = "Raw counts",
  "srs" = "SRS Normalized",
  "css" = "CSS Normalized",

  "optspace" = "Raw+RPCA (2)",
  "gemelli_2" = "RCLR+RPCA (2)",
  "gemelli_3" = "RCLR+RPCA (3)",
  "gemelli_5" = "RCLR+RPCA (5)",
  "gemelli_9" = "RCLR+RPCA (9)",

  "zclr" = "Zero CLR",
  "nclr" = "Nudge CLR",
  "oclr" = "One CLR",
  "deseq2" = "DeSeq2",
  "relabund" = "Relative abundance",
  "sub" = "Subsample",
  "rare" = "Rarefaction",
  "bray_pcoa" = "Bray-Curtis PCoA"
)

pretty_distances <- c(
  euclidean = "Euclidean",
  bray = "Bray-Curtis"
)

pretty_datasets <- c(
  crc = "Colorectal cancer",
  sponges = "Sponges", 
  sleep_apnea = "Sleep Apnea"
)

case_study_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(
    dataset %in% c("crc", "sponges", "sleep_apnea") &
    process != "gemelli_auto" & 
    knn_approach == "knn_pdsgrid" &
    depth != 1000) %>%
    select(-depth, -knn_approach) %>%
  summarize(rep_test = median(test),
            .by = c(dataset, nsamples_per_group, rep, process, distance))
  

pretty_tibble <- case_study_data %>%
  slice_max(nsamples_per_group, by = dataset, with_ties = FALSE) %>%
  select(dataset, nsamples_per_group) %>%
  inner_join(
    tribble(~dataset, ~full_name,
          "crc",     "Colorectal cancer",
          "sponges", "Sponges", 
          "sleep_apnea", "Sleep Apnea"
          ), by = "dataset"
  ) %>%
  mutate(full_name = glue("**{full_name}**<br><span style='font-size:8pt;'>(N={nsamples_per_group} samples per group)</span>"))

case_study_data %>%
  mutate(
    nsamples_per_group = factor(
      nsamples_per_group,
      levels = sort(unique(as.numeric(nsamples_per_group)))),
    dataset = factor(dataset, levels = names(pretty_datasets)),
    process = factor(
      process,
      levels = rev(names(pretty_processes)),
      labels = rev(pretty_processes)
    ),
    distance = factor(
      distance,
      levels = names(pretty_distances))
  ) %>%
  slice_max(nsamples_per_group, by = dataset) %>%
  ggplot(aes(y = process, x = rep_test)) +
  stat_summary(geom = "pointrange", fun.data = median_hilow) +
  facet_grid(distance ~ dataset, scales = "free_y", space = "free_y",
      labeller = labeller(distance = pretty_distances, dataset = pretty_datasets)
  ) +
  labs(
    x = "Accuracy",
    y = NULL, fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black"),
    strip.background = element_rect(color = NA),
    strip.text = element_markdown(size = 12)
  )

ggsave('figures/pointrange_knn_cs.png', width = 8, height = 5)
