#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

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

case_study_data <- read_tsv("data/pooled_data/permanova_treatment.tsv") %>%
  filter(
    dataset %in% names(pretty_datasets) &
    process != "gemelli_auto") %>%
  mutate(
    depth_chr = if_else(depth == 1000, "low", "high")
  ) %>%
  select(-depth, -f_statistic, -r2)

balanced <- case_study_data %>%
  filter(nsamples_per_group != "all") %>%
  summarize(
    f_sig = sum(p_value < 0.05),
    .by = c(dataset, nsamples_per_group, process, distance, depth_chr)
  ) %>%
  filter(depth_chr == "high") %>%
  mutate(
    nsamples_per_group = factor(
      nsamples_per_group,
      levels = unique(as.numeric(nsamples_per_group))),
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
  
balanced_filtered <- balanced %>%
  filter(f_sig != 10)

balanced %>%
  ggplot(aes(x = nsamples_per_group, y = process, fill = f_sig)) +
  geom_tile() +
  geom_text(
    data = balanced_filtered,
    mapping = aes(label = f_sig),
    size = 7, size.unit = "pt") +
  facet_grid(
    distance~dataset, scales = "free", space = "free_y",
    labeller = labeller(
      distance = pretty_distances, dataset = pretty_datasets)
    ) +
  scale_fill_gradient(
    low = "white", high = "dodgerblue",
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
    strip.text.x = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave('figures/heatmap_permanova_sig_cs.png', width = 8, height = 5)


## Move following lines to qmd file? ------
case_study_data %>%
  filter(nsamples_per_group == "all") %>%
  filter(p_value > 0.05)  # RPCA2 and RPCA3 non-significant for CRC
                          # Everything significant for Sponges and Sleep 
                          # Apnea

## the "high" sequence depth produced more significant tests, larger r2 
## and f_statistic values
case_study_data %>%
  filter(nsamples_per_group != "all") %>%
  summarize(
    f_sig = sum(p_value < 0.05),
    .by = c(dataset, nsamples_per_group, process, distance, depth_chr)
  ) %>%
  pivot_wider(names_from = "depth_chr", values_from = "f_sig") %>%
  mutate(diff = high - low) %>%
  filter(diff != 0) %>%
  print(n = Inf)

f <- read_tsv("data/pooled_data/permanova_treatment.tsv") %>%
  filter(
    dataset %in% names(pretty_datasets) &
    process != "gemelli_auto" & 
    nsamples_per_group == "all"
  ) %>%
  mutate(
    depth_chr = if_else(depth == 1000, "low", "high")
  ) %>%
  select(-nsamples_per_group, -rep, -depth, -p_value, -r2) %>% 
  pivot_wider(names_from = "depth_chr", values_from = "f_statistic")
  
f %>%
  ggplot(aes(x = low, y = high, color= dataset)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal() +
  scale_x_log10() +
  scale_y_log10()

ggsave('f_statistic.png')
f %>% filter(low > high)
f %>% filter(process == "rare" & distance ==  "bray")


r2 <- read_tsv("data/pooled_data/permanova_treatment.tsv") %>%
  filter(
    dataset %in% names(pretty_datasets) &
    process != "gemelli_auto" & 
    nsamples_per_group == "all"
  ) %>%
  mutate(
    depth_chr = if_else(depth == 1000, "low", "high")
  ) %>%
  select(-nsamples_per_group, -rep, -depth, -p_value, -f_statistic) %>% 
  pivot_wider(names_from = "depth_chr", values_from = "r2")
  
r2 %>%
  ggplot(aes(x = low, y = high, color= dataset)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal()

ggsave('r2.png')
r2 %>% filter(low > high)
r2 %>% filter(process == "rare" & distance ==  "bray")


case_study_data %>%
  filter(nsamples_per_group == "all" & depth_chr == "high") %>%
  summarize(
    f_sig = sum(p_value < 0.05),
    .by = c(dataset, nsamples_per_group, process, distance, depth_chr)
  ) %>%
  filter(f_sig != 1)

## Move preceding lines to qmd file? ------
