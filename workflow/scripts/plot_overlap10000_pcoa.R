#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

pcoa_fnames <- c("data/simulation/overlap.10000.raw.bray.pcoa.axes",
                 "data/simulation/overlap.10000.oclr.euclidean.pcoa.axes",
                 "data/simulation/overlap.10000.rpca.euclidean.pcoa.axes")

design_fname <- "data/simulation/simulation.design"
output_fname <- "figures/overlap10000.pcoa.png"

plot_pcoa <- function(fname, title, invert_x = FALSE, invert_y = FALSE) {

  pcoa <- read_tsv(fname, col_types = cols(.default = col_double(),
                                           sample = col_character()))

  design <- read_tsv(design_fname, col_types = cols(.default = col_character()))

  loadings <- tibble(percent_exp = colnames(pcoa)[2:3],
                     axis = c(1, 2)) %>%
    mutate(percent_exp = str_replace(percent_exp, ".*_", ""),
          percent_exp = round(as.numeric(percent_exp), digits = 0),
           stylized = glue("Axis {axis} ({percent_exp}%)")) %>%
    pull(stylized)

  ordination <- inner_join(pcoa, design, by = "sample") %>%
    rename_with(~str_replace(., "_.*", ""), contains = "_")

  if (invert_x) {
    ordination <-  mutate(ordination, axis1 = -1 * axis1)
  }
  if (invert_y) {
    ordination <-  mutate(ordination, axis2 = -1 * axis2)
  }


  ordination %>%
    ggplot(aes(x = axis1, y = axis2, color = treatment)) +
    geom_point(size = 1, show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       breaks = c("A", "B"),
                       values = c("tomato", "dodgerblue")) +
    labs(x = loadings[1],
         y = loadings[2],
         title = title) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"))

}

bray <- plot_pcoa(pcoa_fnames[1], "Bray-Curtis")
oclr <- plot_pcoa(pcoa_fnames[2], "Aitchison", invert_x = TRUE)
rpca <- plot_pcoa(pcoa_fnames[3], "Robust PCA")

bray + oclr + rpca +
  plot_layout(widths = 1)

ggsave("figures/overlap10000.pcoa.png", width = 8, height = 3)
