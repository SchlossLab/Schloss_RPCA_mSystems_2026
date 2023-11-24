#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

pcoa_fnames <- c("data/simulation/positive_raw.4000.bray.pcoa.axes",
                 "data/simulation/positive_rpca.4000.euclidean.pcoa.axes",
                 "data/simulation/negative_raw.4000.bray.pcoa.axes",
                 "data/simulation/negative_rpca.4000.euclidean.pcoa.axes")

design_fname <- "data/simulation/overlap.design.gz"
output_fname <- "figures/control4000.pcoa.png"

plot_pcoa <- function(fname, title) {

  pcoa <- read_tsv(fname, col_types = cols(.default = col_double(),
                                           sample = col_character()))

  design <- read_tsv(design_fname, col_types = cols(.default = col_character()))

  loadings <- tibble(percent_exp = colnames(pcoa)[2:3],
                     axis = c(1, 2)) %>%
    mutate(percent_exp = str_replace(percent_exp, ".*_", ""),
           stylized = glue("Axis {axis} ({percent_exp}%)")) %>%
    pull(stylized)

  inner_join(pcoa, design, by = "sample") %>%
    rename_with(~str_replace(., "_.*", ""), contains = "_") %>%
    ggplot(aes(x = axis1, y = axis2, color = treatment)) +
    geom_point(size = 2, show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       breaks = c("A", "B"),
                       values = c("tomato", "dodgerblue")) +
    labs(x = loadings[1],
         y = loadings[2],
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))

}

bray_positive <- plot_pcoa(pcoa_fnames[1], "Bray-Curtis")
rpca_positive <- plot_pcoa(pcoa_fnames[2], "Robust PCA")
bray_negative <- plot_pcoa(pcoa_fnames[3], "Bray-Curtis")
rpca_negative <- plot_pcoa(pcoa_fnames[4], "Robust PCA")

(bray_positive + rpca_positive) / (bray_negative + rpca_negative) +
  plot_layout(widths = 1, heights = 1)

ggsave("figures/control4000.pcoa.png", width = 7, height = 7)