#!/usr/bin/env Rscript

library(tidyverse)
library(ggtext)
library(glue)

study_names <- c("sponges" = "Sponges",
                 "sleep_apnea" = "Sleep Apnea",
                 "crc" = "Colorectal Cancer")

depths <- read_tsv(glue("data/{names(study_names)}/data.depths.tsv"),
                   id = "study",
                   col_types = cols(n_seqs = col_double(),
                                    .default = col_character())) %>%
  mutate(study = str_replace(study, "data/(.*)/.*", "\\1"))

cutoffs <- tibble(study = factor(study_names),
                  published = c(1000, 1000, 0.1),
                  revised = c(13755, 11879, 10439))

treatment_names <- tibble(
  study_treatment = c("crc-healthy", "crc-srn",
                      "sleep_apnea-air", "sleep_apnea-ihh",
                      "sponges-healthy", "sponges-stressed"),
  label = c("Healthy", "SRN",
            "Air", "IHH",
            "Healthy", "Stressed"),
  color = c("A", "B", "A", "B", "A", "B")
)

treatment_names <- depths %>%
  count(study, treatment) %>%
  mutate(study_treatment = glue("{study}-{treatment}")) %>%
  inner_join(treatment_names, by = "study_treatment") %>%
  mutate(label = glue("{label} (N={n})")) %>%
  select(study_treatment, label, color)




depths %>%
  mutate(study_treatment = glue("{study}-{treatment}"),
         study = factor(study,
                        levels = names(study_names),
                        labels = study_names)) %>%
  inner_join(treatment_names, by = "study_treatment") %>%
  select(-sample) %>%
  mutate(label = factor(label, levels = rev(treatment_names$label))) %>%
  ggplot(aes(x = n_seqs, y = label, color = color)) +
  geom_vline(data = cutoffs, aes(xintercept = published), color = "gray",
             linetype = "dashed") +
  geom_vline(data = cutoffs, aes(xintercept = revised), color = "gray") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.9),
              size = 0.75, show.legend = FALSE) +
  facet_wrap(~study, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 1e4, 1e5, 1e6),
                labels = c("10<sup>0</sup>", "10<sup>1</sup>", "10<sup>2</sup>",
                           "10<sup>3</sup>", "10<sup>4</sup>", "10<sup>5</sup>",
                           "10<sup>6</sup>")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_manual(name = NULL,
                     breaks = c("A", "B"),
                     values = c("dodgerblue", "tomato")) +
  coord_cartesian(xlim = c(1, 1e6)) +
  labs(x = "Number of sequences",
       y = NULL) +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 12),
    strip.background = element_blank(),
    axis.text.x = element_markdown(),
    axis.text = element_markdown(color = "black")
  )

ggsave("figures/cs_seq_depth.tiff", width = 5, height = 5,
      compression = "lzw+p")
