#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(broom)
library(cowplot)


read_oclr_data <- function(file_name) {

  condition <- str_split(file_name, pattern = "\\/")[[1]][2]

  fread(file_name, colClasses = c(sample = "character")) |>
    pivot_longer(-sample, names_to = "otu", values_to = "oclr") |>
    mutate(condition = condition)

}


wilcox <- fread("data/overlap_pds/data.10000.100_1.raw.tsv.gz") |>
  pivot_longer(-sample, names_to = "otu", values_to = "count") |>
    mutate(sample = as.numeric(str_replace(sample, "sample_", "")),
          otu = as.numeric(str_replace(otu, "OTU_", "")),
          treatment = if_else(sample <= 100, "A", "B")) %>%
  arrange(sample, otu) %>%
  nest(data = -otu) %>%
  mutate(t = map(.x = data, ~wilcox.test(.x$count ~ .x$treatment) %>% tidy())) %>%
  unnest(t) %>%
  select(otu, p.value) %>% 
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  filter(p.adj < 0.05) %>%
  mutate(condition = "Overlap",
        otu = if_else(otu > 9980, 21 + (otu - 9980), otu))


# skip "otu21" to add a gap in the middle of the heatmap
oclr_data <- c("data/negative_pds/data.10000.100_1.oclr.tsv.gz",
  "data/overlap_pds/data.10000.100_1.oclr.tsv.gz",
  "data/positive_pds/data.10000.100_1.oclr.tsv.gz") |>
    map_dfr(read_oclr_data) |>
    mutate(sample = as.numeric(str_replace(sample, "sample_", "")),
          otu = as.numeric(str_replace(otu, "OTU_", ""))) |>
    filter(otu < 21 | otu > 9980) |>
    mutate(otu = if_else(otu > 9980, 21 + (otu - 9980), otu)) |>
    arrange(condition, sample, otu)

low <- min(oclr_data$oclr)
high <- max(oclr_data$oclr)
mid <- 1 - 2 * high / (high - low)

color_line <- tibble(x = rep(0.25, 200),
                    y = 1:200,
                    treatment = c(rep("A", 100), rep("B", 100)))

heat_map <- oclr_data |>
  mutate(condition = recode(condition, "negative_pds" = "Negative",
                                      "overlap_pds" = "Overlap",
                                      "positive_pds" = "Positive")) |>
  ggplot(aes(x = otu, y = -sample, fill = oclr)) +
  geom_tile() +
  geom_line(data = color_line, inherit.aes = FALSE, show.legend = FALSE,
            aes(x = x, y = -y, color = treatment), linewidth = 3) +
  geom_point(data = wilcox, aes(x = otu), y = -100, inherit.aes = FALSE,
            fill = "white", shape = 23, size = 3) +
  geom_vline(xintercept = 6.5, color = c(NA, "yellow", NA), linewidth = 1) +
  geom_vline(xintercept = 35.5, color = c(NA, "yellow", NA), linewidth = 1) +
  # geom_rect(xmin = 1, xmax = 6, ymin = 20, ymax = , color = "white") +
  scale_fill_gradientn(name = NULL,
                      colors = c("#8E0152", "#FFFFFF", "#276419"),
                      values = scales::rescale(c(-1, mid, 1)),
                      limits = c(low, high)
                      ) +
  scale_color_manual(name = NULL,
                  breaks = c("A", "B"),
                  values = c("tomato", "dodgerblue")) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Features (N=10,000 [40 shown])", y = "Samples (N=200)") +
  facet_grid(condition ~ ., switch = "y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank()
      )

pcoa <- c("data/negative_pds/data.10000.100_1.raw.bray.pcoa.axes",
"data/overlap_pds/data.10000.100_1.raw.bray.pcoa.axes",
"data/positive_pds/data.10000.100_1.raw.bray.pcoa.axes",
"data/negative_pds/data.10000.100_1.gemelli_2.euclidean.pcoa.axes",
"data/overlap_pds/data.10000.100_1.gemelli_2.euclidean.pcoa.axes",
"data/positive_pds/data.10000.100_1.gemelli_2.euclidean.pcoa.axes")


read_pcoa <- function(file_name) {

  simulation <- str_split(file_name, "\\/")[[1]][2]
  metric <- paste(str_split(file_name, "\\.")[[1]][4:5], collapse = "_")

  pcoa_data <- read_tsv(file_name) |>
    mutate(sample = as.numeric(str_replace(sample, "sample_", ""))) |>
    arrange(sample) |>
    mutate(treatment = c(rep("A", 100), rep("B", 100)))

  loadings <- colnames(pcoa_data)[2:3] |>
    str_replace("axis[12]_", "") |>
    as.numeric()

  colnames(pcoa_data) <- c("sample", "axis1", "axis2", "treatment")

  tibble(axis1 = loadings[1], axis2 = loadings[2],
        simulation = simulation, metric = metric,
        nest(pcoa_data))

}

pcoa_combined <- map_dfr(pcoa, read_pcoa) |>
  mutate(simulation = recode(simulation, "negative_pds" = "Negative",
                                         "overlap_pds" = "Overlap",
                                         "positive_pds" = "Positive"),
        metric = recode(metric, "raw_bray" = "Bray-Curtis",
                                "gemelli_2_euclidean" = "RPCA"))

pcoa_data <- pcoa_combined |>
  select(simulation, metric, data) |>
  unnest(data)

loadings <- pcoa_combined |>
  select(simulation, metric, axis1, axis2) |>
  mutate(x_min = min(pcoa_data$axis1), x_max = max(pcoa_data$axis1),
         y_min = min(pcoa_data$axis2), y_max = max(pcoa_data$axis2),
         axis1 = format(round(axis1, 1), 1),
         axis2 = format(round(axis2, 1), 1))


pcoa_figure <- pcoa_data |>
  ggplot(aes(x = axis1, y = axis2, color = treatment)) +
  geom_point(show.legend = FALSE) +
  geom_text(data = loadings, aes(x = x_min, y = y_max, label = axis2),
            inherit.aes = FALSE, size = 2, hjust = 0.25) +
  geom_text(data = loadings, aes(x = x_max, y = y_min, label = axis1),
            inherit.aes = FALSE, size = 2, hjust = 0.75) +
  facet_grid(metric ~ simulation) +
  coord_fixed() +
  labs(x = "PCo Axis 1", y = "PCo Axis 2") +
  scale_color_manual(name = NULL,
                    breaks = c("A", "B"),
                    values = c("tomato", "dodgerblue")) +
  theme_classic() +
  theme(
    panel.background = element_rect(color = "black", linewidth = 1.1),
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black")
  )


combined <- plot_grid(heat_map, pcoa_figure,
                      labels = "AUTO", label_size = 18, hjust = 0,
                      ncol = 1, rel_heights = c(0.65, 0.85))

ggsave(filename = "figures/heatmap_pcoa_pds.tiff",
      combined,
      width = 7, height = 5.75,
      compression = "lzw+p")
