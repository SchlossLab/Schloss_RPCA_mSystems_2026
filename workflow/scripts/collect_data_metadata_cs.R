#!/usr/bin/env Rscript

library(tidyverse)
library(biomformat)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

# biom_fname <- "data/sponges/table.biom"
# metadata_fname <- "data/sponges/metadata.tsv"
# column <- "health_status"

biom_fname <- args[1]
metadata_fname <- args[2]
column <- args[3]

output_path <- tolower(
  str_replace(biom_fname,
              "data/(.*)/table.biom",
              "\\1")
)

count_data <- read_hdf5_biom(biom_fname) %>%
  biom() %>%
  .[["data"]] %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  mutate(otu = glue("OTU_{1:nrow(.)}"), .before = everything()) %>%
  pivot_longer(cols = - otu, names_to = "sample", values_to = "count")

metadata <- read_tsv(metadata_fname) %>%
  select(sample = sample_name, treatment = all_of(column)) %>%
  mutate(treatment = tolower(treatment))


inner_join(count_data, metadata, by = "sample") %>%
  pivot_wider(id_cols = c(sample, treatment),
              names_from = otu, values_from = count) %>%
  write_tsv(glue("data/{output_path}/data_treatment.tsv.gz"))
