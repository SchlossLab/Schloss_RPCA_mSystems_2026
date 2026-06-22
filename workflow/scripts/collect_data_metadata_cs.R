#!/usr/bin/env Rscript

library(tidyverse)
library(biomformat)
library(glue)
library(ape) 

args <- commandArgs(trailingOnly = TRUE)
biom_fname <- args[1]
metadata_fname <- args[2]
column <- args[3]

# biom_fname <- "data/sponges/table.biom"
# metadata_fname <- "data/sponges/metadata.tsv"
# column <- "health_status"

output_path <- tolower(
  str_replace(biom_fname,
              "data/(.*)/table.biom",
              "\\1")
)

biom_data <- read_hdf5_biom(biom_fname) %>%
  biom()

otu_names <- biom_data[["rows"]] %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  unnest(id) %>%
  pull(id)

count_data <- biom_data %>%
  .[["data"]] %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  mutate(otu = otu_names, .before = everything()) %>%
  pivot_longer(cols = - otu, names_to = "sample", values_to = "count")


metadata <- read_tsv(metadata_fname) %>%
  select(sample = sample_name, treatment = all_of(column)) %>%
  mutate(treatment = tolower(treatment))

count_metadata <- inner_join(count_data, metadata, by = "sample") %>%
  mutate(n = sum(count), .by = "otu") %>%
  filter(n != 0) %>% 
  pivot_wider(id_cols = c(sample, treatment),
              names_from = otu, values_from = count)
              
count_metadata %>%
  write_tsv(glue("data/{output_path}/data_treatment.tsv.gz"))

otu_names <- colnames(count_metadata)[-c(1, 2)]

glue(">{otu_names}\n{otu_names}") %>%
  write_lines(glue("data/{output_path}/data_treatment.fasta"))
