# simulate data for CENSORING CASE
do_one_sim <- function(n_obs = 1000,        # sample size
                       iter,
                       sim_type = c("simple", "hvtn505"),
                       effect_type = c("moderate", "strong", "null"),
                       path = "/global/scratch/nhejazi/txshift-sims/data/"
                      ) {
  # default to simple simulation setting
  sim_type <- match.arg(sim_type)
  effect_type <- match.arg(effect_type)

  if (sim_type == "simple") {
    # baseline covariate -- simple, binary
    W <- cbind(rnorm(n_obs, 3, sd = 1), rbinom(n_obs, 1, 0.6),
               rbinom(n_obs, 1, 0.6 / 2))

    ## create treatment based on binary baseline W
    A <- rnorm(n_obs, mean = 2 * rowSums(W[, -1]), sd = 1)

    # create outcome as a function of A, W + white noise
    Y <- rbinom(n_obs, 1, as.numeric(plogis(rowMeans(W) - A)))

    # censoring (multi-stage sample) process
    C <- rbinom(n_obs, 1, plogis(W[, 2] + Y))

    # subset to observed data
    X <- data.table::as.data.table(cbind(W, A, Y, C))
    data.table::setnames(X, c(paste0("W", seq_len(ncol(W))), "A", "Y", "Delta"))

    # save full data
    saveRDS(X, file = paste0(path,
                             "sim_nsamp_", n_obs,
                             "_dataset_", iter, ".rds"))

  } else if (sim_type == "hvtn505") {
    ## baseline covariates, similar to those in HVTN 505
    W <- cbind(rnorm(n_obs, mean = 26.6, sd = 5.7),  # mimics BMI
               rpois(n_obs, lambda = 40),            # mimics age
               rbinom(n_obs, 1, 0.4),                # mimics behavioral risk
               rbinom(n_obs, 1, 0.3))                # mimics race (categorical)

    ## create treatment based on binary baseline W
    A <- -1.37 + 0.015 * W[, 2] + 0.004 * W[, 1] + 0.05 * W[, 3] +
      0.25 * W[, 4] + rnorm(n_obs, mean = 0, sd = 0.2)

    ## create outcome as a function of A, W + white noise
    if (effect_type == "moderate") {
      Y <- rbinom(n_obs, 1, as.numeric(plogis(-2.8 - 0.0016 * W[, 2] +
                                              0.0013 * W[, 1] + 0.0678 *
                                              W[, 3] + 0.039 * W[, 4] -
                                              0.033 * A)))
    } else if (effect_type == "strong") {
      Y <- rbinom(n_obs, 1, as.numeric(plogis(-4.8 - 0.0016 * W[, 2] +
                                              0.0013 * W[, 1] + 0.0678 *
                                              W[, 3] + 0.039 * W[, 4] -
                                              3.3 * A)))
    } else {
      Y <- rbinom(n_obs, 1, as.numeric(plogis(-2.9 - 0.0016 * W[, 2] +
                                              0.0013 * W[, 1] + 0.0678 *
                                              W[, 3] + 0.039 * W[, 4])))
    }

    # two-phased sampling process (includes all cases)
    C <- rbinom(n_obs, 1, plogis(-2.45 - 0.027 * W[, 1]  + 0.166 * W[, 4] +
                                 0.012 * W[, 2] + 0.39 * W[, 3]))
    C[Y == 1] <- 1

    # subset to observed data
    X <- data.table::as.data.table(cbind(W, A, Y, C))
    data.table::setnames(X, c(paste0("W", seq_len(ncol(W))), "A", "Y", "Delta"))

    # save full data
    saveRDS(X, file = paste0(path,
                             "sim_nsamp_", n_obs,
                             "_effect_", effect_type,
                             "_dataset_hvtn505_", iter, ".rds"))
  }
}
