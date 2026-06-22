#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)

simulations <- c("overlap_pds")

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

overlap_power <- read_tsv("data/pooled_data/permanova_treatment.tsv.gz") %>%
  filter(dataset %in% simulations) %>%
  filter(process %in% processes) %>%
  select(-dataset, -nsamples_per_group) %>%
  group_by(depth, process, distance) %>%
  summarize(power = mean(p_value < 0.05), .groups = "drop") %>%
  select(depth, process, distance, metric = power) %>%
  mutate(analysis = "power") %>%
  rbind(tibble(depth = c(1000, 2000), process = c("zclr", "zclr"),
              distance = c("euclidean", "euclidean"), metric = c(-1, -1),
              analysis = c("power", "power")))

overlap_knn <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(dataset %in% simulations) %>%
  filter(process %in% processes) %>%
  filter(knn_approach == "knn_pdsgrid") %>%
  select(-dataset, -nsamples_per_group, -knn_approach, -k_size) %>%
  group_by(depth, process, distance, rep) %>%
  summarize(metric = median(test), .groups = "drop") %>%
  mutate(analysis = "knn")

overlap_knn <- overlap_knn %>%
  rbind(tibble(depth = c(1000, 2000), process = c("zclr", "zclr"),
              distance = c("euclidean", "euclidean"), rep = c(1,1),
              metric = c(-1, -1),
              analysis = c("knn", "knn")))

vline <- tibble(distance = c("euclidean", "bray"),
                analysis = c("knn", "knn"),
                xintercept = 0.5)

overlap_plot <- bind_rows(overlap_knn, overlap_power) %>%
  mutate(depth = fct_reorder(as.character(depth), -depth),
        distance = factor(distance,
                          levels = c("euclidean", "bray"),
                          labels = c("Euclidean", "Bray-Curtis")),
        process = factor(process, levels = processes),
        analysis = factor(analysis, levels = c("power", "knn"),
                           labels = c("Power",
                                      "Accuracy"))) %>%
  ggplot(aes(y = process, x = metric, fill = depth)) +
  geom_vline(xintercept = 0.5, color = rep(c("white", "gray"), 2)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.95),
              position = position_dodge(width = 0.75), shape = 21,
              orientation = "y") +
  facet_grid(distance ~ analysis,
            scales = "free", space = "free_y", switch = "x") +
  coord_cartesian(xlim = c(0, 1)) +
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
    strip.text.y = element_text(size = 10, color = "black"),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = NA, color = "black",
                                    linewidth = 1.1),
    legend.text = element_text(hjust = 1),
    legend.key.height = unit(0.8, "lines"),
    legend.key = element_rect(color = NA),
    panel.ontop = TRUE,
    panel.spacing = unit(10, "pt")
  )

labeled_overlap_plot <- ggdraw(overlap_plot) +
  draw_label("PERMANOVA", x = 0.34, y = 0.97,
            fontface = "bold", size = 12) +
  draw_label("KNN Classification", x = 0.635, y = 0.97,
            fontface = "bold", size = 12)

ggsave(filename = "figures/fst_knn_pds.tiff",
      labeled_overlap_plot,
      width = 6, height = 5,
      compression = "lzw+p")
