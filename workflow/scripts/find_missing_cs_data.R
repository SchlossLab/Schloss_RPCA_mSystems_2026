#!/usr/bin/env Rscript

library(tidyverse)
library(biomformat)
library(glue)
library(ape) 



process_data <- function(dataset, column){

  biom_fname <- glue("data/{dataset}/table.biom")
  tree_fname <- glue("deicode_benchmarking/case_studies/data/{dataset}/tree_relabelled.tre.gz")
  metadata_fname <- glue("data/{dataset}/metadata.tsv")

  # read in biom formatted file and connect with representative sequence; make
  # count data tidy
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


  # read in metadata file and get appropriate treatment column
  metadata <- read_tsv(metadata_fname) %>%
    select(sample = sample_name, treatment = all_of(column)) %>%
    mutate(treatment = tolower(treatment))


  # connect count data with metadata; remove representative sequences that no
  # longer are in dataset
  count_metadata <- inner_join(count_data, metadata, by = "sample") %>%
    mutate(n = sum(count), .by = "otu") %>%
    filter(n != 0) %>% 
    pivot_wider(id_cols = c(sample, treatment),
                names_from = otu, values_from = count)
                

  # read in tree and remove anything that isn't in the otu data
  tree <- read.tree(tree_fname)
  otu_names <- colnames(count_metadata)[-c(1, 2)]
  remove_names <- tree$tip.label[!(tree$tip.label %in% otu_names)]

  tree <- tree %>%
    drop.tip(remove_names)


  # find what in otu data is not found in tree
  missing_names_in_tree <- otu_names[!(otu_names %in% tree$tip.label)]
  n_missing_names_in_tree <- length(missing_names_in_tree)
  n_otus_total <- length(otu_names)
  n_tree_total <- length(tree$tip.label)

  percent_missing <- count_metadata %>%
    select(-treatment) %>%
    pivot_longer(-sample) %>%
    mutate(missing = name %in% missing_names_in_tree) %>%
    summarize(removed = sum(value), .by = missing) %>%
    mutate(percent_removed = removed / sum(removed) * 100) %>%
    filter(missing) %>%
    pull(percent_removed)

  tibble(
    n_otus = n_otus_total,
    n_tree_tips = n_tree_total,
    n_missing = n_missing_names_in_tree,
    percent_missing = percent_missing)
}

tibble(
  dataset = c("sponges", "sleep_apnea"),
  column = c("health_status", "exposure_type")
) %>%
  mutate(analysis = map2(dataset, column, process_data)) %>%
  unnest(analysis) %>%
  write_tsv("data/pooled_data/missing_cs_data.tsv")

# from this code analyzing data from before removing poorly sequenced samples
#   dataset     column        n_otus n_tree_tips n_missing percent_missing
#   <chr>       <chr>          <int>       <int>     <int>           <dbl>
# 1 sponges     health_status  11160       11074        86           0.689
# 2 sleep_apnea exposure_type   1369        1345        24          39.7  


# from code analyzing data from after removing poorly sequenced samples
#   dataset     n_otus n_tree_tips n_missing percent_missing
#   <chr>        <int>       <int>     <int>           <dbl>
# 1 sponges      11101       11017        84           0.686
# 2 sleep_apnea   1157        1137        20          39.8  