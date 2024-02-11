#!/usr/bin/env Rscript

# this code was inspired largely from
# workflow/scripts/collect_data_metadata_sim.py and makes use of
# workflow/scripts/simulations.R. I have increased the number of features to
# 10000 from 1000 to allow for a power analysis that gives a moderate level of
# power with partial overlap of features across ranks.


library(tidyverse)
source("workflow/scripts/simulations.R")

args <- commandArgs(trailingOnly = TRUE)
desired_depth <- as.integer(args[1])
rep <- as.integer(args[2])

# desired_depth <- 4000
# rep <- 1

set.seed(rep)

n_samples <- 200
n_features <- 10000
n_ranks <- 2
nsamples_tag <- as.integer(n_samples / 2)

# the number to the right of the comment was the original version in the ipynb
# and script indicated below. I selected larger values that would reliably give
# a minimum sequencing depth per sample at or above the desired sequencing
# depth. I won't do anything with 500 since it wasn't used in the paper and give
# a 0 negative control matrix (would need to add 1 to the mean abundance)

seq_depth <- c("500" = 3.05e2,
               "1000" = 750, #6.1e2,
               "2000" = 1500, #1.25e3,
               "4000" = 2800, #2.5e3,
               "10000" = 7250) #6.05e3

fitting_parameter <- seq_depth[[as.character(desired_depth)]]

################################################################################

generate_counts_data <- function(depth, parameter, scaling_factor, overlap,
                                nrank, nsamples, notus) {

  raw_counts <- block_diagonal_gaus_features(nsamples, notus, nrank,
                                    overlap, parameter) %>%
    homoscedastic(parameter / scaling_factor) %>%
    heteroscedastic(parameter / scaling_factor) %>%
    subsample(parameter)

  min_sample_sum <- min(apply(raw_counts, 2, sum))

  stopifnot(min_sample_sum > depth)

  raw_counts %>%
    t() %>%
    as_tibble(.name_repair = \(x) glue("OTU_{1:ncol(.)}")) %>%
    mutate(sample = glue("sample_{1:nrow(.)}"), .before = everything())

}

################################################################################
# now we want to subsample each sample to desired_depth number of sequences and
# make the data frame wide again with samples in rows and otus in columns
################################################################################

subsample_counts_data <- function(mat, depth){

  mat %>%
    pivot_longer(-sample) %>%
    filter(value > 0) %>%
    uncount(value) %>%
    slice_sample(n = depth, by = sample) %>%
    count(sample, name) %>%
    pivot_wider(names_from = name, values_from = n, values_fill = 0)

}

################################################################################
#
# Positive control - subsampled to desired depth
#
################################################################################

overlap <- 0

generate_counts_data(desired_depth, fitting_parameter, 60, overlap,
                                n_ranks, n_samples, n_features) %>%
  subsample_counts_data(desired_depth) %>%
  write_tsv(glue("data/positive_pds/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz"))


################################################################################
#
# Negative control - basically random noise, subsampled to desired depth
#
################################################################################

matrix(0, nrow = n_samples, ncol = n_features) %>%
    homoscedastic(fitting_parameter / 60) %>%
    heteroscedastic(fitting_parameter / 60) %>%
    subsample(fitting_parameter) %>%
    as_tibble(.name_repair = \(x) glue("OTU_{1:ncol(.)}")) %>%
    mutate(sample = glue("sample_{1:nrow(.)}"), .before = everything()) %>%
    subsample_counts_data(desired_depth) %>%
    write_tsv(glue("data/negative_pds/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz"))


################################################################################
#
# Overlap control
#
################################################################################

overlap <- 4994

generate_counts_data(desired_depth, fitting_parameter, 15, overlap,
                                n_ranks, n_samples, n_features) %>%
  subsample_counts_data(desired_depth) %>%
  write_tsv(glue("data/overlap_pds/data.{desired_depth}.{nsamples_tag}_{rep}.raw.tsv.gz"))

################################################################################
#
# Design files
#
################################################################################

design <- tibble(
  sample = glue("sample_{1:n_samples}"),
  treatment = rep(c("A", "B"), each = nsamples_tag)
)

write_tsv(
  design,
  glue("data/positive_pds/data.{desired_depth}.{nsamples_tag}_{rep}.design")
)

write_tsv(
  design,
  glue("data/negative_pds/data.{desired_depth}.{nsamples_tag}_{rep}.design")
)

write_tsv(
  design,
  glue("data/overlap_pds/data.{desired_depth}.{nsamples_tag}_{rep}.design")
)
