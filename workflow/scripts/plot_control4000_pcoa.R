#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

pcoa_fnames <- c("data/positive/data.4000.100_1.raw.bray.pcoa.axes",
                 "data/positive/data.4000.100_1.gemelli_2.euclidean.pcoa.axes",
                 "data/negative/data.4000.100_1.raw.bray.pcoa.axes",
                 "data/negative/data.4000.100_1.gemelli_2.euclidean.pcoa.axes")

design_fnames <- c("data/positive/data.4000.100_1.design",
                  "data/negative/data.4000.100_1.design")

output_fname <- "figures/control4000.pcoa.png"

plot_pcoa <- function(axes_fname, design_fname, title) {

  pcoa <- read_tsv(axes_fname, col_types = cols(.default = col_double(),
                                           sample = col_character()))

  design <- read_tsv(design_fname, col_types = cols(.default = col_character()))

  loadings <- tibble(percent_exp = colnames(pcoa)[2:3],
                     axis = c(1, 2)) %>%
    mutate(percent_exp = str_replace(percent_exp, ".*_", ""),
          percent_exp = round(as.numeric(percent_exp), digits = 0),
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
         subtitle = title) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(face = "bold"),
          plot.title.position = "plot"
          )

}

bray_positive <- plot_pcoa(pcoa_fnames[1], design_fnames[1], "Bray-Curtis") +
  labs(title = "Positive control")
rpca_positive <- plot_pcoa(pcoa_fnames[2], design_fnames[1], "Robust PCA")

bray_negative <- plot_pcoa(pcoa_fnames[3], design_fnames[2], "Bray-Curtis") +
  labs(title = "Negative control")
rpca_negative <- plot_pcoa(pcoa_fnames[4], design_fnames[2], "Robust PCA")

(bray_positive + rpca_positive) / (bray_negative + rpca_negative) +
  plot_layout(widths = 1, heights = 1)

ggsave(output_fname, width = 7, height = 7)
