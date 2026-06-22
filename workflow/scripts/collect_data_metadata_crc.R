#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ape)

shared_fname <- "data/crc/data.otu.shared.gz"
metadata_fname <- "data/crc/sra_info.tsv.gz"
tree_fname <- "data/crc/data.otu.tre"
fasta_fname <- "data/crc/data.otu.fasta"

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


################################################################################
###
### Need to relabel tree and count_table with OTU rather than sequence number
###
################################################################################

seq_to_otu_lookup <- readLines(fasta_fname) %>%
  str_subset(">") %>% 
  str_remove(">") %>%
  str_remove("\\|.*") %>%
  as_tibble_col() %>%
  separate_wider_delim(value, delim = "\t", names = c("sequence", "otu")) %>%
  mutate(otu = Otu_to_OTU(otu))


tree <- read.tree(tree_fname)
tree$tip.label <- seq_to_otu_lookup[["otu"]][match(tree$tip.label,
                                               seq_to_otu_lookup[["sequence"]])]
write.tree(tree, "data/crc/data_treatment.tre")
