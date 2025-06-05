#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(cowplot)
library(RColorBrewer)

simulations <- c("overlap_pds")

processes <- c("raw", "optspace", "gemelli_2", "zclr", "nclr", "oclr", "deseq2",
               "css", "srs", "bray_pcoa")

pretty_processes <- c(
  "raw" = "Raw",
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

overlap_knn <- read_tsv("data/pooled_data/knn_results.tsv.gz") %>%
  filter(dataset %in% simulations) %>%
  filter(process %in% processes) %>%
  select(-dataset, -nsamples_per_group, -k_size) %>%
  group_by(knn_approach, depth, process, distance, rep) %>%
  summarize(metric = median(test), .groups = "drop")

overlap_knn <- overlap_knn %>%
  rbind(tibble(depth = rep(c(1000, 2000), 3),
              process = rep(c("zclr", "zclr"), 3),
              distance = rep(c("euclidean", "euclidean"), 3),
              rep = rep(c(1,1), 3),
              metric = rep(c(-1, -1), 3),
              knn_approach = rep(c("knn_cmfixed", "knn_cmgrid", "knn_pdsgrid"), each = 2)))

overlap_plot <- overlap_knn %>%
  mutate(depth = fct_reorder(as.character(depth), -depth),
        distance = factor(distance,
                          levels = c("euclidean", "bray"),
                          labels = c("Euclidean", "Bray-Curtis")),
        process = factor(process, levels = processes),
        knn_approach = factor(
                        knn_approach,
                        levels = c("knn_cmfixed", "knn_cmgrid", "knn_pdsgrid"),
                        labels = c("60:40 split\nw/ K=13",
                                   "60:40 split\nw/ Optimized K",
                                   "80:20 split\nw/ Optimized K")
                      )) %>%
  ggplot(aes(y = process, x = metric, fill = depth)) +
  geom_vline(xintercept = 0.5, color = "gray") +
  stat_summary(fun = median, fun.min = min, fun.max = max,
              position = position_dodge(width = 0.75), shape = 21,
              orientation = "y") +
  facet_grid(distance ~ knn_approach,
            scales = "free", space = "free_y") +
  coord_cartesian(xlim = c(0.4, 1)) +
  scale_fill_manual(breaks = c(1000, 2000, 4000, 10000),
                    values = RColorBrewer::brewer.pal(name = "OrRd", n = 4),
                    labels = c("1,000", "2,000", "4,000", "10,000")) +
  scale_y_discrete(breaks = processes, labels = pretty_processes[processes],
                  limits = rev) +
  labs(x = "Accuracy", y = NULL, fill = NULL,
      ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 10),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    panel.spacing = unit(10, "pt"),
    panel.background = element_rect(fill = NA, color = "black",
                                    linewidth = 1.1),
    legend.text = element_text(hjust = 1),
    legend.key = element_rect(color = NA),
    legend.key.height = unit(0.8, "lines"),
    panel.ontop = TRUE
  )

ggsave(filename = "figures/knn_methods_pds.tiff",
      overlap_plot,
      width = 7, height = 5,
      compression = "lzw+p")
