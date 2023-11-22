#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

pcoa_fnames <- c("data/simulation/overlap_raw.10000.bray.pcoa.axes",
                 "data/simulation/overlap_rpca.10000.euclidean.pcoa.axes")

design_fname <- "data/simulation/overlap_raw.design.gz"
output_fname <- "figures/overlap10000.pcoa.png"

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
                       values = c("red", "dodgerblue")) +
    labs(x = loadings[1],
         y = loadings[2],
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))

}

bray <- plot_pcoa(pcoa_fnames[1], "Bray-Curtis + PCoA")
rpca <- plot_pcoa(pcoa_fnames[2], "Robust PCA")

bray + rpca +
  plot_layout(widths = 1)

ggsave("figures/overlap10000.pcoa.png", width = 7, height = 3.5)