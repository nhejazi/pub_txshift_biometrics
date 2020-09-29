get_truth <- function(n_obs = 1e7,     # sample size
                      w_mean = 3,      # mean of W when Gaussian
                      w_prob = 0.6,    # proportion of W=1 in group
                      tx_mult = 2,     # multiplier effect of W=1 on treatment
                      delta) {
  # baseline covariate -- simple, binary
  W <- cbind(rnorm(n_obs, w_mean, sd = 1), rbinom(n_obs, 1, w_prob),
             rbinom(n_obs, 1, w_prob / 2))

  ## create treatment based on binary baseline W
  A <- rnorm(n_obs, mean = tx_mult * rowSums(W[, -1]), sd = 1)

  # create outcome as a function of A, W + white noise
  Qbar_Aplusdelta_W <- as.numeric(plogis(rowMeans(W) - (A + delta)))
  return(mean(Qbar_Aplusdelta_W))
}

get_truth_hvtn505 <- function(n_obs = 1e7, delta,
                              type = c("moderate", "null", "strong")) {
  ## baseline covariates, similar to those in HVTN 505
  W <- cbind(rnorm(n_obs, mean = 26.6, sd = 5.7),  # mimics BMI
             rpois(n_obs, lambda = 40),            # mimics age
             rbinom(n_obs, 1, 0.4),                # mimics behavioral risk
             rbinom(n_obs, 1, 0.3))                # mimics race (categorical)

  ## create treatment based on binary baseline W
  A <- -1.37 + 0.015 * W[, 2] + 0.004 * W[, 1] + 0.05 * W[, 3] +
    0.25 * W[, 4] + rnorm(n_obs, mean = 0, sd = 0.2)

  ## create outcome as a function of A, W + white noise
  if (type == "moderate") {
    Qbar_Aplusdelta_W <- rbinom(n_obs, 1,
                                as.numeric(plogis(-2.8 - 0.0016 * W[, 2] +
                                           0.0013 * W[, 1] + 0.0678 *
                                           W[, 3] + 0.039 * W[, 4] -
                                           0.033 * (A + delta))))
  } else if (type == "null") {
    Qbar_Aplusdelta_W <- rbinom(n_obs, 1,
                                as.numeric(plogis(-2.9 - 0.0016 * W[, 2] +
                                           0.0013 * W[, 1] + 0.0678 *
                                           W[, 3] + 0.039 * W[, 4])))
  } else if (type == "strong") {
    Qbar_Aplusdelta_W <- rbinom(n_obs, 1,
                                as.numeric(plogis(-4.8 - 0.0016 * W[, 2] +
                                           0.0013 * W[, 1] + 0.0678 *
                                           W[, 3] + 0.039 * W[, 4] -
                                           3.3 * (A + delta))))
  }
  return(mean(Qbar_Aplusdelta_W))
}
