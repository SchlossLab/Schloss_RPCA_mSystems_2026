#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)

simulations <- c("overlap")

processes <- c("raw", "optspace", "gemelli_2", "zclr", "nclr", "oclr", "deseq2",
               "css", "srs", "bray_pcoa")

pretty_processes <- c(
  "raw" = "Raw counts",
  "srs" = "SRS Normalized",
  "css" = "CSS Normalized",

  "optspace" = "Raw + RPCA",
  "gemelli_2" = "RCLR + RPCA",
  "zclr" = "Zero CLR",
  "nclr" = "Nudge CLR",
  "oclr" = "One CLR",
  "deseq2" = "DeSeq2",
  "bray_pcoa" = "Bray-Curtis PCoA"
)

pretty_distance <- c(
  euclidean = "Euclidean",
  bray = "Bray-Curtis"
)

overlap_fst <- read_tsv("data/pooled_data/permanova_treatment.tsv.gz") %>%
  filter(dataset %in% simulations) %>%
  filter(process %in% processes) %>%
  select(-dataset, -nsamples_per_group) %>%
  select(depth, process, distance, rep, metric = f_statistic) %>%
  mutate(analysis = "f_st")

overlap_knn <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(dataset %in% simulations) %>%
  filter(process %in% processes) %>%
  filter(knn_approach == "knn_pdsgrid") %>%
  select(-dataset, -nsamples_per_group, -knn_approach, -k_size) %>%
  group_by(depth, process, distance, rep) %>%
  summarize(metric = median(test), .groups = "drop") %>%
  mutate(analysis = "knn")

vline <- tibble(distance = c("euclidean", "bray"),
                analysis = c("knn", "knn"),
                xintercept = 0.5)

overlap_plot <- bind_rows(overlap_knn, overlap_fst) %>%
  mutate(depth = fct_reorder(as.character(depth), -depth),
        distance = factor(distance,
                          levels = c("euclidean", "bray"),
                          labels = c("Euclidean", "Bray-Curtis")),
        process = factor(process, levels = processes),
        analysis = factor(analysis, levels = c("f_st", "knn"),
                           labels = c("F-Statistic",
                                      "Accuracy"))) %>%
  ggplot(aes(y = process, x = metric, fill = depth)) +
  geom_vline(xintercept = 0.5, color = rep(c("white", "gray"), 2)) +
  stat_summary(fun = median, fun.min = min, fun.max = max,
              position = position_dodge(width = 0.75), shape = 21,
              orientation = "y") +
  facet_grid(distance ~ analysis,
            scales = "free", space = "free_y", switch = "x") +
  coord_cartesian(xlim = c(0, NA)) +
  scale_fill_manual(breaks = c(1000, 2000, 4000, 10000),
                    values = RColorBrewer::brewer.pal(name = "OrRd", n = 4),
                    labels = c("1,000", "2,000", "4,000", "10,000")) +
  scale_y_discrete(breaks = processes, labels = pretty_processes[processes],
                  limits = rev) +
  labs(x = NULL, y = NULL, fill = NULL,
      subtitle = "" #Opens up some space for draw_label
      ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 10),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = NA, color = "black",
                                    linewidth = 1.1),
    legend.text = element_text(hjust = 1),
    legend.key = element_rect(color = NA),
    legend.key.height = unit(0.8, "lines"),
    panel.ontop = TRUE
  )

labeled_overlap_plot <- ggdraw(overlap_plot) +
  draw_label("PERMANOVA", x = 0.34, y = 0.97,
            fontface = "bold", size = 12) +
  draw_label("KNN Classification", x = 0.635, y = 0.97,
            fontface = "bold", size = 12)

ggsave(filename = "figures/fst_knn_martino.tiff",
      labeled_overlap_plot,
      width = 6, height = 5,
      compression = "lzw+p")
