#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

shared_fname <- "data/crc/data.otu.shared.gz"
metadata_fname <- "data/crc/sra_info.tsv.gz"

Otu_to_OTU <- function(x){
  str_replace(x, "Otu0*", "OTU_")
}

count_data <- fread(shared_fname,
                    colClasses = list(character = "Group")) %>%
  select(sample = Group, starts_with("Otu")) %>%
  rename_with(Otu_to_OTU, starts_with("Otu"))

metadata <- read_tsv(metadata_fname,
                     col_types = cols(Sample_Name = col_character())) %>% 
  mutate(treatment = if_else(diagnosis %in% c("Cancer", "adv Adenoma"),
                             "srn",
                             "healthy")) %>%
  select(sample = Sample_Name, treatment) %>%
  distinct()

inner_join(count_data, metadata, by = "sample") %>%
  select(sample, treatment, starts_with("OTU")) %>%
  write_tsv("data/crc/data_treatment.tsv.gz")
