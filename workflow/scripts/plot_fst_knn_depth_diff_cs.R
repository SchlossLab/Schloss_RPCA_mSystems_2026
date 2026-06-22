#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(ggtext)

permanova <- read_tsv("data/pooled_data/permanova_treatment.tsv.gz")

n_seqs_cs <- permanova %>%
  filter(nsamples_per_group == "all" & depth != 1000 & process == "raw" & distance == "bray") %>%
  select(dataset, depth) %>%
  mutate(depth = format(depth, big.mark = ",")) %>%
  pivot_wider(names_from = dataset, values_from = depth)


pretty_processes <- c(
  "raw" = "Raw counts",
  "css" = "CSS Normalized",
  "srs" = "SRS Normalized",

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
  bray = "Bray-Curtis",
  wunifrac = "Weighted\nUniFrac"
)

# Need to insert the Ns
pretty_datasets <- c(
  sponges = glue("Sponges\n(N={n_seqs_cs$sponges})"), 
  sleep_apnea = glue("Sleep Apnea\n(N={n_seqs_cs$sleep_apnea})"),
  crc = glue("Colorectal Cancer\n(N={n_seqs_cs$crc})")
)

pretty_metric <- c(
  fst = "**PERMANOVA**<br>(F statistic)",
  knn = "**KNN Classification**<br>(Accuracy)"
)

knn_data <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(knn_approach == "knn_pdsgrid",
         dataset %in% names(pretty_datasets),
         process != "gemelli_auto") %>%
  slice_max(nsamples_per_group, by = dataset) %>%
  mutate(depth_chr = if_else(depth == 1000, "low", "high")) %>%
  select(dataset, depth_chr, process, distance,
        rep, iteration, test) %>%
  summarize(test = median(test), .by = c(dataset, depth_chr, process, distance, rep)) %>%
  pivot_wider(names_from = "depth_chr", values_from = "test") %>%
  mutate(diff = 100 * (high - low) / low) 

fst_data <- permanova %>%
  filter(dataset %in% names(pretty_datasets),
         nsamples_per_group == "all",
         process != "gemelli_auto") %>%
  mutate(depth_chr = if_else(depth == 1000, "low", "high")) %>%
  select(dataset, depth_chr, process, distance, f_statistic) %>%
  pivot_wider(names_from = "depth_chr", values_from = "f_statistic") %>%
  mutate(diff = 100 * (high - low) / low, rep = 1) 


bind_rows(knn = knn_data, fst = fst_data, .id = "metric") %>%
  filter(
    !(dataset == "sleep_apnea" & process == "bray_pcoa" & metric != "knn")
  ) %>%
  mutate(
    # dataset = factor(
    #   dataset,
    #   levels = names(pretty_datasets),
    #   labels = pretty_datasets),
    process = factor(
      process,
      levels = rev(names(pretty_processes)),
      labels = rev(pretty_processes)),
    distance = factor(
      distance,
      levels = names(pretty_distances),
      labels = pretty_distances),
    metric = factor(
      metric,
      levels = names(pretty_metric),
      labels = pretty_metric)
  ) %>%
  ggplot(
    aes(x = diff, y = process, color = dataset)) +
    geom_vline(xintercept = 0, color = "gray", linewidth = 0.3) +
    stat_summary(fun = median, fun.min = min, fun.max = max,
                position = position_dodge(width = 0.5), size = 0.3) + 
    facet_grid(distance~metric, scales = "free", space = "free_y") +
    labs(
      x = "Increase in parameter when using N\nsequences relative to 1,000 sequences (%)",
      y = NULL,
      title = NULL,
      subtitle = NULL,
      caption = NULL
    ) +
    scale_color_manual(
      name = NULL,
      breaks = names(pretty_datasets),
      labels = pretty_datasets,
      values = c('#1b9e77','#d95f02','#7570b3')
    ) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = NA, color = "black",
                                     linewidth = 1.1),
      panel.ontop = TRUE,

      axis.text = element_text(color = "black"),
      axis.title.x = element_text(margin = margin(t = 10)),

      strip.background = element_blank(),
      strip.text.x = element_markdown(size = 12),
      strip.text.y = element_text(size = 10),
      strip.clip = "off",

      legend.key = element_blank(),
      legend.box.spacing = unit(-5, "pt"),
      legend.key.spacing.y = unit(5, "pt"),
      legend.text = element_text(size = 8)
    )

ggsave(
  "figures/fst_knn_depth_diff_cs.tiff", width = 7, height = 7,
  compression = "lzw+p")
