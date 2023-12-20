#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

depth <- args[1]
seed <- args[2]

# depth <- 10000
# seed <- 8

dir.create("figures/pcoa", showWarnings = FALSE)

output_fname <- glue("figures/pcoa/{depth}_{seed}.tiff")

process_distance <- c("optspace.euclidean" = "OptSpace",
                      "raw.euclidean" = "Euclidean",
                      "raw.bray" = "Bray-Curtis",
                      "oclr.euclidean" = "Aitchison",
                      "gemelli_2.euclidean" = "RPCA")

simulations <- c("negative" = "Negative",
                "overlap" = "Overlap",
                "positive" = "Positive")


plot_pcoa <- function(fname) {

  design_fname <- str_replace(fname, "(.*100_\\d*)\\..*", "\\1.design")

  pcoa <- read_tsv(fname,
                   col_types = cols(.default = col_double(),
                                    sample = col_character()))

  design <- read_tsv(design_fname,
                     col_types = cols(.default = col_character()))

  loadings <- tibble(percent_exp = colnames(pcoa)[2:3],
                     axis = c(1, 2)) %>%
    mutate(percent_exp = str_replace(percent_exp, ".*_", ""),
          percent_exp = round(as.numeric(percent_exp), digits = 0),
           stylized = glue("Axis {axis} ({percent_exp}%)")) %>%
    pull(stylized)

  ordination <- inner_join(pcoa, design, by = "sample") %>%
    rename_with(~str_replace(., "_.*", ""), contains = "_")


  simulation <- str_replace(fname, "data/(.*)/data.*", "\\1")
  distance <- str_replace(fname,
                          "data/.*/data\\.\\d*\\.100_\\d*\\.(.*).pcoa.axes",
                          "\\1")

  title <- if_else(distance == "optspace.euclidean",
                   simulations[simulation],
                   NA_character_)

  tag <- if_else(simulation == "negative",
                 process_distance[distance],
                 NA_character_)

  p <- ordination %>%
    ggplot(aes(x = axis1, y = axis2, color = treatment)) +
    geom_point(size = 1, show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       breaks = c("A", "B"),
                       values = c("dodgerblue", "tomato")) +
    labs(x = loadings[1],
         y = loadings[2]) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 10)
    )

  if (!is.na(title)) {
    p <- p + labs(title = title)
  }

  if (!is.na(tag) && !is.na(title)) {
    p <- p +
        labs(tag = tag) +
        theme(plot.tag = element_text(angle = 90, hjust = 0.5,
                                      face = "bold",
                                      margin = margin(b = -145, l = 0, r = 10)
                                      )
              )
  } else if (!is.na(tag)) {
    p <- p +
        labs(tag = tag) +
        theme(plot.tag = element_text(angle = 90, hjust = 0.5,
                                      face = "bold",
                                      margin = margin(b = -105, l = 0, r = 10)
                                      )
              )
  }

  return(p)
}


grid_plots <- expand_grid(sim = names(simulations),
                          proc_dist = names(process_distance)) %>%
  mutate(
    sim = factor(sim, levels = names(simulations)),
    proc_dist = factor(proc_dist, levels = names(process_distance)),
    path = glue("data/{sim}/data.{depth}.100_{seed}.{proc_dist}.pcoa.axes")
  ) %>%
  mutate(plots = map(path, plot_pcoa)) %>%
  arrange(proc_dist, sim)

(grid_plots$plots[[1]]  + grid_plots$plots[[2]]  + grid_plots$plots[[3]])  /
(grid_plots$plots[[4]]  + grid_plots$plots[[5]]  + grid_plots$plots[[6]])  /
(grid_plots$plots[[7]]  + grid_plots$plots[[8]]  + grid_plots$plots[[9]])  /
(grid_plots$plots[[10]] + grid_plots$plots[[11]] + grid_plots$plots[[12]]) /
(grid_plots$plots[[13]] + grid_plots$plots[[14]] + grid_plots$plots[[15]]) +
  plot_layout(widths = 1, heights = 1)

ggsave(output_fname,
      width = 7.5, height = 10,
      compression = "lzw+p")
