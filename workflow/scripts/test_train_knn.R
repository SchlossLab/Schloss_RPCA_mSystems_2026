library(tidyverse)
library(glue)


partition_samples <- function(m, frac = 0.8, label) {

  train_label <- glue("train_{label}")

  m %>%
    nest(data = -treatment) %>%
    mutate("{train_label}" := map(data,
                            function(x) {
                              n_train <- as.integer(frac * nrow(x))
                              sample(
                                c(rep("train", n_train),
                                  rep(label, nrow(x) - n_train))
                              )
                            }
                          )) %>%
  select(treatment, all_of(train_label), data) %>%
  unnest(c(all_of(train_label), data))

}

get_validation_accuacy <- function(m, k, truth, v_fraction) {

  m %>%
    partition_samples(frac = v_fraction, label = "validate") %>%
    get_accuracy(k, "validate", truth)

}

get_optimal_ksize <- function(m, k_grid, truth, n_splits, v_fraction) {

  expand_grid(replicate = 1:n_splits,
              k_size = k_grid) %>%
    mutate(accuracy = map_dfr(.data$k_size,
                              ~get_validation_accuacy(m, .x, truth, v_fraction))) %>%
    unnest(accuracy) %>%
    group_by(k_size) %>%
    summarize(validate = mean(validate)) %>%
    slice_max(order_by = validate, with_ties = FALSE)

}

get_accuracy <- function(m, k, validate_test, truth) {

  train_column <- paste0("train_", validate_test)

  m %>%
    filter(get(train_column) == "train") %>%
    select(-starts_with("train_"), -all_of(.[["sample"]])) %>%
    pivot_longer(cols = -c(treatment, sample), names_to = validate_test) %>%
    rename(train = sample) %>%
    group_by_at(validate_test) %>%
    slice_min(order_by = value, n = k, with_ties = FALSE) %>%
    summarize(prop_case = sum(treatment) / n()) %>%
    mutate(classification = prop_case > 0.5) %>%
    inner_join(., truth, by = setNames("sample", validate_test)) %>%
    summarize("{validate_test}" := sum(classification == treatment) / n())

}

run_one_knn_loop <- function(m, k_grid, v_fraction, t_fraction, n_splits) {

  truth <- m %>%
    select(sample, treatment)

  train_test <- m %>%
    partition_samples(frac = t_fraction, label = "test")

  train <- train_test %>%
    filter(train_test == "train") %>%
    select(treatment, sample, all_of(.[["sample"]]))

  get_optimal_ksize(train, k_grid, truth, n_splits, v_fraction) %>%
    mutate(test = get_accuracy(train_test, k = k_size, "test", truth)) %>%
    unnest(test)

}

run_knn <- function(design_fname, dmatrix_fname, output_fname,
                    v_fraction = 0.8, t_fraction = 0.8,
                    n_iloops = 100, n_oloops = 100,
                    k_grid = seq(1, 100, 10)) {

  set.seed(1)

  k_grid <- as.integer(k_grid)

  treatments <- read_tsv(design_fname,
                        col_types = cols(.default = col_character()))

  case <- treatments[[1, "treatment"]]

  dmatrix <- read_tsv(dmatrix_fname,
                      col_types = cols(sample = col_character(),
                                      .default = col_double()))

  treat_dmatrix <- inner_join(treatments, dmatrix, by = "sample") %>%
                  mutate(treatment = treatment == case)

  tibble(iteration = 1:n_oloops) %>%
    mutate(knn_loops = map_dfr(.data$iteration,
                              ~run_one_knn_loop(treat_dmatrix, k_grid,
                                                v_fraction, t_fraction,
                                                n_iloops))) %>%
    unnest(knn_loops) %>%
    write_tsv(output_fname)

}
