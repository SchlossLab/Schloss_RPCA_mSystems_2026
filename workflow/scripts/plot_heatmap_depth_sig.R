library(tidyverse)

# NOTE: We don't really care about the treatment effect in this simulation. It's
# p-value will be inflated since we're artificially duplicating every sample.
# This will make each treatment look more similar to itself than it really is,
# which will reduce the intra-treatment variation increasing the F statistic
# and decreasing the p-value for the treatment effect.

pretty_processes <- c(
  "raw" = "Raw counts",
  "srs" = "SRS Normalized",
  "css" = "CSS Normalized",

  "optspace" = "Raw+RPCA (2)",
  "gemelli_2" = "RCLR+RPCA (2)",
  "gemelli_3" = "RCLR+RPCA (3)",
  "gemelli_5" = "RCLR+RPCA (5)",
  "gemelli_9" = "RCLR+RPCA (9)",

  "zclr" = "Zero CLR",
  "nclr" = "Nudge CLR",
  "oclr" = "One CLR",
  "deseq2" = "DeSeq2",
  "relabund" = "Relative abundance",
  "sub" = "Subsample",
  "rare" = "Rarefaction",
  "bray_pcoa" = "Bray-Curtis PCoA"
)

pretty_distances <- c(
  euclidean = "Euclidean",
  bray = "Bray-Curtis"
)

pretty_datasets <- c(
  crc = "CRC",
  sponges = "Sponges", 
  sleep_apnea = "Sleep\nApnea",
  negative = "Negative",
  overlap = "Overlap",
  positive = "Positive"
)

pretty_experiments <- c(
  original = "Original\nsimulations",
  revised = "Revised\nsimulations",
  case_study = "Case\nstudies"
)


read_tsv("data/pooled_data/permanova_size.tsv") %>%
#  filter(rep == 1) %>%
  filter(effect == "size" & process != "gemelli_auto") %>% 
  select(dataset, rep, process, distance, p_value) %>%
  summarize(n_sig = sum(p_value < 0.05), .by = c(dataset, process, distance)) %>%
  mutate(
    experiment = case_when(dataset %in% c("negative", "positive", "overlap") ~ "original",
                           dataset %in% c("crc", "sleep_apnea", "sponges") ~ "case_study",
                           TRUE ~ "revised"),
    experiment = factor(experiment,
                      levels = names(pretty_experiments),
                      labels = pretty_experiments),
    distance = factor(distance,
                      levels = names(pretty_distances),
                      labels = pretty_distances),
    dataset = str_replace(dataset, "_pds", ""),
    dataset = factor(dataset,
                      levels = names(pretty_datasets),
                      labels = pretty_datasets),
    process = factor(process,
                    levels = rev(names(pretty_processes)),
                    labels = rev(pretty_processes))
  ) %>%
  ggplot(aes(y = process, x = dataset, fill = n_sig)) +
  geom_tile() +
  facet_grid(distance~ experiment, scales = "free", space = "free_y") +
  scale_fill_gradient(
    low = "white", high = "dodgerblue",
    n.breaks = 6
  ) +

  labs(x = NULL, y = NULL, fill = NULL) +
  coord_cartesian(expand = FALSE) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(color = NA),
    strip.text = element_text(size = 12),
    strip.text.x = element_text(face = "bold")
  )


ggsave('figures/heatmap_depth_sig.png', width = 8, height = 5)
