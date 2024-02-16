#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
set.seed(19760620)

# This function is an R version of the original python code for the
# block_diagonal_gaus function found in the original benchmarking dataset
# (deicode_benchmarking/simulations/scripts/simulations.py). The problems with
# this function include: (i) between the two ranks, the samples overlap rather
# than the features and (ii) the abundance distribution for each feature varies
# by the sample rather than by the feature
block_diagonal_gaus_samples <- function(ncols, nrows, nblocks, overlap,
                                       minval, maxval) {

    mat <- matrix(0, nrow = nrows, ncol = ncols)

    block_cols <- as.integer(ncols / nblocks)
    block_rows <- as.integer(nrows / nblocks)

    gradient <- 5
    mu <- seq(0, 10, length.out = block_cols + overlap)
    sigma <- 2.0

    b <- matrix(rep(dnorm(gradient, mu, sigma), block_rows),
                nrow = block_rows, ncol = block_cols + overlap,
                byrow = TRUE) * maxval

    mat[1:block_rows, 1:(block_cols + overlap)] <- b

    # Make last block, fill in the remainder
    gradient <- 5
    mu <- seq(0, 10, length.out = block_cols + overlap)
    sigma <- 4

    b <- matrix(rep(dnorm(gradient, mu, sigma), block_rows),
                nrow = block_rows, ncol = block_cols + overlap,
                byrow = TRUE) * maxval

    mat[(block_rows + 1):nrows, (block_cols - overlap + 1):ncols] <- b

    return(mat)
}

# The following function is a revised version of the previous function that
# overcomes the two problems with the original. Namely, here the features
# overlap between the ranks and the abundance distribution of the features
# varies by feature rather than by sample within the same rank. I have largely
# kept the spirit and parameterization of the original function. One notable
# difference is that the abundance of the overlapping features draws from a
# common distribution with a mean zero and a moderately sized standard
# deviation; any resulting values less than zero are set to zero.

block_diagonal_gaus_features <- function(ncols, nrows, nblocks, overlap,
                                       maxval) {

    mat <- matrix(0, nrow = nrows, ncol = ncols)

    block_cols <- as.integer(ncols / nblocks)
    block_rows <- as.integer(nrows / nblocks)

    gradient <- 5
    mu <- seq(0, 10, length.out = block_rows - overlap)
    sigma <- 4 #match second block

    b <- matrix(rep(dnorm(gradient, mu, sigma), block_cols),
                nrow = block_rows - overlap, ncol = block_cols,
                byrow = FALSE) * maxval

    mat[1:(block_rows - overlap), 1:block_cols] <- b

    # Make last block, fill in the remainder
    gradient <- 5
    mu <- seq(0, 10, length.out = block_rows - overlap)
    sigma <- 4

    b <- matrix(rep(dnorm(gradient, mu, sigma), block_cols),
                nrow = block_rows - overlap, ncol = block_cols,
                byrow = FALSE) * maxval

    mat[(block_rows + overlap + 1):nrows, (block_cols + 1):ncols] <- b

    if(overlap != 0) {
    # This line will fill the overlapping features with random numbers using
    # a mean of 0 and standard deviation of 50. If the resulting value is less
    # than 0, then the value is returned as zero. This helps add some sparisty.
    # The value of the mean could be increased or decreased to affect the level
    # of sparsity in the overlapping region.
      rand_value <- rnorm(2 * overlap * ncols, mean = 0, sd = 50)
      mat[(block_rows - overlap + 1):(block_rows + overlap), ] <- 
        if_else(rand_value < 0, 0, rand_value)
    }

    return(mat)
}



homoscedastic <- function(mat, intensity) {

    nr <- nrow(mat)
    nc <- ncol(mat)

    matrix(rnorm(length(mat), mean = mat, sd = intensity),
          nrow = nr, ncol = nc)

}

heteroscedastic <- function(mat, intensity) {

    i <- sample(size = 5000, nrow(mat), replace = TRUE)
    j <- sample(size = 5000, ncol(mat), replace = TRUE)

    for (index in 1:5000) {
      mat[i[index], j[index]] <- rnorm(1,
                                      mean = mat[i[index], j[index]],
                                      sd = intensity)
    }

    mat[mat < 0] <- 0 #set all negative values to zero

    return(mat)

}

subsample <- function(mat, sparsity) {

    #""" yij ~ PLN( lambda_{ij}, /phi ) """
    mu <- sparsity * t(t(mat) / apply(mat, 2, sum))

    #in following if mu value is zero, then the result is zero
    matrix(rpois(n = length(mat), rlnorm(n = length(mat), meanlog = log(mu))),
            nrow = nrow(mat), ncol = ncol(mat))

}
