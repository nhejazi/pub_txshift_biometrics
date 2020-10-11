# use custom package library
.libPaths("/global/scratch/nhejazi/R")

# packages and options
library(here)
library(foreach)
library(future)
library(doFuture)
library(doRNG)
library(data.table)
library(tidyverse)
library(hal9001)
library(haldensify)
library(sl3)
library(txshift)

# parallelization
options(future.globals.maxSize = 10^24)
registerDoFuture()
plan(multiprocess)
set.seed(43972)

# parameters to specify simulation scenario
n_sim <- 2000
estim <- c("tmle", "onestep")
sim_type <- "simple"
#sim_type <- "hvtn505"

# parameters that vary across simulation scenarios
if (sim_type == "simple") {
  shift_grid <- seq(-0.5, 0.5, by = 0.5)
  n_samp <- cumsum(rep(sqrt(100), 5))^2
  fit_type_gQ <- "mle"
  param_type <- "shift"
  #param_type <- "npmsm"
} else if (sim_type == "hvtn505") {
  #shift_grid <- seq(-2, 2, by = 0.5)
  #n_samp <- 1400
  #fit_type_gQ <- "hal"
  #param_type <- "shift"
  #param_type <- "npmsm"
  #effect_type <- "moderate"
  #effect_type <- "strong"
  #effect_type <- "null"
}

# load scripts to set up simulation and estimation procedure
source(here("R", "utils", "compute_estimators.R"))

# loop over estimator types for convenience
for (estim_type in estim) {
  # assess all sample sizes by looping over
  sim_results <- map(n_samp, function(sample_size) {
    # simulation for a given sample size
    result_nsamp <- foreach(this_iter = seq_len(n_sim),
                            .options.multicore = list(preschedule = FALSE),
                            .errorhandling = "remove") %dorng% {
      sim_iter_out <- estim_one_sim(n_obs = sample_size,
                                    delta_shift = shift_grid,
                                    pi_fit = "glm",
                                    ipcw_eff = TRUE,
                                    estimator = estim_type,
                                    iter = this_iter,
                                    sim_type = sim_type,
                                    effect_type = effect_type,
                                    param_type = param_type,
                                    fit_type_gQ = fit_type_gQ)
      if (param_type == "shift") {
        sim_iter_out <- sim_iter_out %>%
          add_column(sim_id = rep(this_iter, length(shift_grid)))
      } else if (param_type == "npmsm") {
        sim_iter_out <- sim_iter_out %>%
          add_column(sim_id = rep(this_iter, 1))
      }
      sim_iter_out
    }
    result_nsamp <- result_nsamp %>%
      bind_rows() %>%
      as_tibble()
    if (param_type == "shift") {
      colnames(result_nsamp) <-
        c("delta", "lwr_ci", "param_est", "upr_ci", "param_var", "sim_id")
    } else if (param_type == "npmsm") {
      colnames(result_nsamp) <-
        c("lwr_ci", "param_est", "upr_ci", "param_se", "p_value", "sim_id")
    }

    # save results to file
    timestamp <- str_replace_all(Sys.time(), " ", "_")
    saveRDS(object = result_nsamp,
            file = here("data", ifelse(sim_type == "simple", "simple",
                                       paste(sim_type, effect_type,
                                             sep = "_")),
                        paste0("02-eff_", estim_type, "-pi_glm_", param_type,
                               "_n", sample_size, "_", timestamp, ".rds"))
           )
  })
}
