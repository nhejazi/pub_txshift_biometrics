# use custom package library
.libPaths("/global/scratch/nhejazi/R")

# packages and options
library(here)
library(future)
library(future.apply)
library(data.table)
library(tidyverse)
options(scipen = 999)
plan(multiprocess, workers = 20)
seed_int <- 43972

# load scripts to set up simulation and estimation procedure
source(here("R", "utils", "make_censored_data.R"))

# set simulation settings
n_sim <- 2000
sim_type <- "simple"
#sim_type <- "hvtn505"
## NOTE: ignored for sim_type != "hvtn505"
effect_type <- "moderate"
#effect_type <- "strong"
#effect_type <- "null"

# setting-dependent simulation parameters
if (sim_type == "simple") {
  n_samp <- (cumsum(rep(sqrt(100), 7))^2)
} else if (sim_type == "hvtn505") {
  n_samp <- 1400
}

# assess all sample sizes by looping over
sim_data <- map(n_samp, function(sample_size) {
  # simulation for a given sample size
  data_nsamp <-
    future.apply::future_lapply(X = seq_len(n_sim), FUN = function(x_iter) {
      sim_iter_out <- do_one_sim(n_obs = sample_size,
                                 iter = x_iter,
                                 sim_type = sim_type,
                                 effect_type = effect_type)
    },
    future.seed = seed_int
  )
})
