#!/usr/bin/env Rscript

library(tidyverse)

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
    dataset %in% names(pretty_datasets) &
    process != "gemelli_auto" & 
    knn_approach == "knn_pdsgrid" &
    depth != 1000) %>%
    select(-depth, -knn_approach) %>%
  summarize(rep_test = median(test),
            .by = c(dataset, nsamples_per_group, rep, process, distance)) %>%
  summarize(test = median(rep_test),
            .by = c(dataset, nsamples_per_group, process, distance)) %>%
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
  )
  
# label <- case_study_data %>%
#   slice_max(nsamples_per_group, by = dataset) %>%
#   mutate(label = format(round(test, 2L), nsmall = 2L))


case_study_data %>%
  mutate(label = format(round(test, 2L), nsmall = 2L)) %>%
  ggplot(aes(x = nsamples_per_group, y = process, fill = test, label = label)) +
  geom_tile(show.legend = FALSE) +
  geom_text(size = 6, size.unit = "pt") +
  facet_grid(
    distance~dataset, scales = "free", space = "free_y",
    labeller = labeller(
      distance = pretty_distances, dataset = pretty_datasets)
    ) +
  scale_fill_gradient(
    low = "white", high = "red",
    n.breaks = 6
  ) +
  labs(
    x = "Number of samples per group",
    y = NULL, fill = NULL
  ) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(size = 12),
    strip.text.x = element_text(face = "bold")
  )

ggsave('figures/heatmap_knn_cs.png', width = 8, height = 5)
