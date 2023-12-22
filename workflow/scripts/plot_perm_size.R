#!/usr/bin/env Rscript

library(tidyverse)
library(glue)

median_iqr <- function(x) {

    tibble(ymin = quantile(x, prob = 0.25),
            y = median(x),
            ymax = quantile(x, prob = 0.75))
}

datasets <- c("negative" = "Negative",
              "overlap" = "Overlap",
              "positive" = "Positive",
              "sponges" = "Sponges",
              "sleep_apnea" = "Sleep apnea",
              "crc" = "CRC")

processes <- c("raw", "rare", "sub", "oclr", "nclr", "zclr", "gemelli_2", "optspace")

read_tsv("data/pooled_data/permanova_size.tsv") %>%
  filter(effect == "size",
        process %in% processes) %>%
  select(dataset, process, process, distance, f_statistic, p_value) %>%
  mutate(process_dist = glue("{process}_{distance}"),
        dataset = factor(dataset, levels = names(datasets), labels = datasets)) %>%
  group_by(dataset, process_dist) %>%
  summarize(f_sig = mean(p_value < 0.05)) %>%
  ggplot(aes(x = dataset, y = process_dist, fill = f_sig)) +
  geom_tile(color = "black", show.legend = FALSE) +
  scale_fill_gradient(name = NULL,
                      low = "#FFFFFF",
                      high = "#FF0000") +
  labs(x = NULL,
      y = NULL) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank()
  )


ggsave("figures/permanova_size.tiff",
      width = 7, height = 3,
      compression = "lzw+p")
