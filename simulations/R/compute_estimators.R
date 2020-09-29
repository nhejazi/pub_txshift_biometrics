# compute estimators from censored/sampled data
estim_one_sim <- function(n_obs = 1000,             # sample size
                          delta_shift = 0.5,        # shift parameter
                          pi_fit = "hal",           # GLM or HAL for TMLE?
                          ipcw_eff = TRUE,          # efficient IPCW-TMLE?
                          estimator =               # TML or one-step estimator
                            c("tmle", "onestep"),
                          naive = FALSE,            # ignore sampling mechanism
                          iter,
                          sim_type = c("simple", "hvtn505"),
                          effect_type = c(TRUE, FALSE),
                          param_type = c("shift", "npmsm"),
                          fit_type_gQ = c("mle", "hal"),
                          path = "/global/scratch/nhejazi/txshift-sims/data/"
                         ) {
  # read in full data
  if (sim_type == "simple") {
    X <- readRDS(file = paste0(path, "sim_nsamp_", n_obs, "_dataset_",
                               iter, ".rds"))
  } else if (sim_type == "hvtn505") {
    X <- readRDS(file = paste0(path, "sim_nsamp_", n_obs, "_effect_",
                               effect_type, "_dataset_", "hvtn505_",
                               iter, ".rds"))
  }
  w_names <- str_subset(colnames(X), "W")

  # data computing sampling mechanism
  V <- cbind(X[, ..w_names], X[, "Y"])

  # build observed data from full data
  O <- X %>%
    dplyr::filter(Delta == 1) %>%
    dplyr::select(-Delta)

  # fit sampling mechanism using HAL
  if (pi_fit == "hal") {
    pi_mech <- fit_hal(Y = X$Delta, X = as.matrix(V),
                       standardize = FALSE,
                       family = "binomial",
                       fit_type = "glmnet",
                       lambda.min.ratio = 1e-5,
                       yolo = FALSE)
    pi_mech_pred <- predict(pi_mech, new_data = V)
  } else {
    subsetX <- c(w_names, "Y", "Delta")
    pi_mech <- glm(Delta ~ ., data = X[, ..subsetX],
                   family = binomial())
    pi_mech_pred <- fitted(pi_mech)
  }
  ipc_weights_out <-
    unname((as.numeric(X$Delta == 1) / pi_mech_pred)[X$Delta == 1])
  ipcw_out <- list(pi_mech = pi_mech_pred, ipc_weights = ipc_weights_out)
  ipcw_fit_args_in <- list(fit_type = "fit_spec")
  ipcw_fit_spec_in <- ipcw_out

  if (fit_type_gQ == "mle") {
    # compute the MLEs for the treatment mechanism densities
    g_data <- O %>%
      dplyr::select(-Y)
    glm_fit <- glm(A ~ ., data = g_data, weights = ipc_weights_out)

    # MLE of beta
    # beta_n <- matrix(glm_fit$coefficients, ncol = 1)
    w_times_beta_n <- fitted(glm_fit)

    # get the estimate of sigma
    ## 1) normalized weights
    ipw_trunc_norm <- ipc_weights_out / sum(ipc_weights_out)

    ## 2) weighted mean of A
    abar <- weighted.mean(O$A, w = ipw_trunc_norm)

    ## 3) weighted total variance of A
    tot_var <- sum((O$A - abar)^2 * ipw_trunc_norm)

    ## 4) variance of conditional mean
    var_cond_mean <- sum((w_times_beta_n -
                          weighted.mean(w_times_beta_n, w = ipw_trunc_norm))^2 *
                         ipw_trunc_norm)

    ## 5) conditional variance of Y|X -> get MLE of sigma
    cond_var <- tot_var - var_cond_mean
    sigma_n <- sqrt(cond_var)

    results_shifted_means <- lapply(delta_shift, function(delta) {
      # compute the densities for the treatment mechanism via the MLE
      gn_downshift <- dnorm(O$A - delta, mean = w_times_beta_n, sd = sigma_n)
      gn_upupshift <- dnorm(O$A + 2 * delta, mean = w_times_beta_n,
                            sd = sigma_n)
      gn_upshift <- dnorm(O$A + delta, mean = w_times_beta_n, sd = sigma_n)
      gn_noshift <- dnorm(O$A, mean = w_times_beta_n, sd = sigma_n)
      gn_mle <- as.data.table(cbind(gn_downshift, gn_noshift, gn_upshift,
                                    gn_upupshift))
      setnames(gn_mle, c("downshift", "noshift", "upshift", "upupshift"))

      # use wrapper function to fit either the TMLE or the one-step estimator
      if (!naive) {
        est_shift <- txshift(W = X[, ..w_names],
                             A = X$A,
                             Y = X$Y,
                             delta = delta,
                             C = X$Delta,
                             V = V,
                             estimator = estimator,
                             max_iter = 5,
                             fluctuation = "weighted",
                             g_fit_args = list(fit_type = "fit_spec"),
                             gn_fit_spec = gn_mle,
                             ipcw_fit_args = ipcw_fit_args_in,
                             ipcw_fit_spec = ipcw_fit_spec_in,
                             Q_fit_args = list(fit_type = "glm",
                                               glm_formula = paste0("Y ~ A + ",
                                                  paste0(w_names,
                                                         collapse = " + "))),
                             ipcw_efficiency = ipcw_eff,
                             eif_reg_type = "hal"              # HAL by default
                            )
      } else {
        O <- data.table::as.data.table(O)
        est_shift <- txshift(W = O[, ..w_names],
                             A = O$A,
                             Y = O$Y,
                             delta = delta,
                             estimator = estimator,
                             max_iter = 5,
                             fluctuation = "weighted",
                             g_fit_args = list(fit_type = "fit_spec"),
                             gn_fit_spec = gn_mle,
                             Q_fit_args = list(fit_type = "glm",
                                               glm_formula = paste0("Y ~ A + ",
                                                  paste0(w_names,
                                                         collapse = " + ")))
                            )
      }
      if (param_type == "shift") {
        # compute inference and create output object
        ci <- confint(est_shift)
        out <- list(delta = delta,
                    lwr_ci = ci[1],
                    est = ci[2],
                    upr_ci = ci[3],
                    var = est_shift$var) %>%
          as_tibble()
        return(out)
      } else if (param_type == "npmsm") {
        return(est_shift)
      }
    })
    if (param_type == "npmsm") {
      # multiplier for CI construction
      ci_level <- 0.95
      ci_mult <- c(1, -1) * stats::qnorm((1 - ci_level) / 2)

      # matrix of EIF(O_i) values and estimates across each parameter estimated
      psi_vec <- sapply(results_shifted_means, `[[`, "psi")
      eif_mat <- sapply(results_shifted_means, `[[`, "eif")

      # set weights to be the inverse of the variance of each TML estimate
      #wts <- as.numeric(1 / diag(stats::cov(eif_mat)))
      wts <- rep(1, length(delta_shift))

      # compute the MSM parameters
      intercept <- rep(1, length(delta_shift))
      x_mat <- cbind(intercept, delta_shift)
      omega <- diag(wts)
      s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
      msm_param <- as.vector(s_mat %*% psi_vec)

      # compute inference for MSM using individual EIF(O_i) for each parameter
      msm_eif <- tcrossprod(eif_mat, s_mat)
      msm_var <- diag(stats::cov(msm_eif))
      msm_se <- sqrt(msm_var / nrow(msm_eif))

      # build confidence intervals and hypothesis tests for EIF(msm)
      ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
      pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

      # create summary table for MSM estimates
      msm_out <- list(
        param = names(msm_se),
        ci_low = ci_msm_param[, 1],
        param_est = msm_param,
        ci_high = ci_msm_param[, 2],
        param_se = msm_se,
        p_value = pval_msm_param
      ) %>%
      as_tibble()
    }
  } else if (fit_type_gQ == "hal") {
    # coerce to data.table
    O <- data.table::as.data.table(O)

    # fit HAL for g and Q outside of loop involving shift parameter
    gn_hal_fit <- haldensify(W = as.matrix(O[, ..w_names]),
                             A = O$A, wts = ipc_weights_out)

    Qn_hal_fit <- fit_hal(Y = O$Y,
                          X = as.matrix(cbind(O[, ..w_names], O$A, O$Y)),
                          standardize = FALSE, family = "binomial",
                          fit_type = "glmnet", lambda.min.ratio = 1e-5,
                          weights = ipc_weights_out, yolo = FALSE)

    results_shifted_means <- lapply(delta_shift, function(delta) {
      # predict from conditional treatment density model for given shift
      gn_hal_downshift <- predict(gn_hal_fit,
                                  new_W = as.matrix(O[, ..w_names]),
                                  new_A = (O$A - delta))
      gn_hal_upupshift <- predict(gn_hal_fit,
                                  new_W = as.matrix(O[, ..w_names]),
                                  new_A = (O$A + 2 * delta))
      gn_hal_upshift <- predict(gn_hal_fit,
                                new_W = as.matrix(O[, ..w_names]),
                                new_A = (O$A + delta))
      gn_hal_noshift <- predict(gn_hal_fit,
                                new_W = as.matrix(O[, ..w_names]),
                                new_A = O$A)
      gn_halmle <- as.data.table(cbind(gn_hal_downshift, gn_hal_noshift,
                                       gn_hal_upshift, gn_hal_upupshift))
      setnames(gn_halmle, c("downshift", "noshift", "upshift", "upupshift"))

      # predict from outcome model
      Qn_noshift <- predict(Qn_hal_fit,
                            new_data = as.matrix(cbind(O[, ..w_names],
                                                       O$A, O$Y)))
      Qn_upshift <- predict(Qn_hal_fit,
                            new_data = as.matrix(cbind(O[, ..w_names],
                                                       O$A + delta, O$Y)))
      Qn_halmle <- as.data.table(cbind(Qn_noshift, Qn_upshift))
      setnames(Qn_halmle, c("noshift", "upshift"))

      # use wrapper function to fit either the TMLE or the one-step estimator
      if (!naive) {
        est_shift <- txshift(W = X[, ..w_names],
                             A = X$A,
                             Y = X$Y,
                             delta = delta,
                             C = X$Delta,
                             V = V,
                             estimator = estimator,
                             max_iter = 5,
                             fluctuation = "weighted",
                             g_fit_args = list(fit_type = "fit_spec"),
                             gn_fit_spec = gn_halmle,
                             ipcw_fit_args = ipcw_fit_args_in,
                             ipcw_fit_spec = ipcw_fit_spec_in,
                             Q_fit_args = list(fit_type = "fit_spec"),
                             Qn_fit_spec = Qn_halmle,
                             ipcw_efficiency = ipcw_eff,
                             eif_reg_type = "hal"              # HAL by default
                            )
      } else {
        O <- data.table::as.data.table(O)
        est_shift <- txshift(W = O[, ..w_names],
                             A = O$A,
                             Y = O$Y,
                             delta = delta,
                             estimator = estimator,
                             max_iter = 5,
                             fluctuation = "standard",
                             g_fit_args = list(fit_type = "fit_spec"),
                             gn_fit_spec = gn_halmle,
                             Q_fit_args = list(fit_type = "fit_spec"),
                             Qn_fit_spec = Qn_halmle
                            )
      }
      if (param_type == "shift") {
        # compute inference and create output object
        ci <- confint(est_shift)
        out <- list(delta = delta,
                    lwr_ci = ci[1],
                    est = ci[2],
                    upr_ci = ci[3],
                    var = est_shift$var) %>%
          as_tibble()
        return(out)
      } else if (param_type == "npmsm") {
        return(est_shift)
      }
    })
    if (param_type == "npmsm") {
      # multiplier for CI construction
      ci_level <- 0.95
      ci_mult <- c(1, -1) * stats::qnorm((1 - ci_level) / 2)

      # matrix of EIF(O_i) values and estimates across each parameter estimated
      psi_vec <- sapply(results_shifted_means, `[[`, "psi")
      eif_mat <- sapply(results_shifted_means, `[[`, "eif")

      # set weights to be the inverse of the variance of each TML estimate
      #wts <- as.numeric(1 / diag(stats::cov(eif_mat)))
      wts <- rep(1, length(delta_shift))

      # compute the MSM parameters
      intercept <- rep(1, length(delta_shift))
      x_mat <- cbind(intercept, delta_shift)
      omega <- diag(wts)
      s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
      msm_param <- as.vector(s_mat %*% psi_vec)

      # compute inference for MSM using individual EIF(O_i) for each parameter
      msm_eif <- tcrossprod(eif_mat, s_mat)
      msm_var <- diag(stats::cov(msm_eif))
      msm_se <- sqrt(msm_var / nrow(msm_eif))

      # build confidence intervals and hypothesis tests for EIF(msm)
      ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
      pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

      # create summary table for MSM estimates
      msm_out <- list(
        param = names(msm_se),
        ci_low = ci_msm_param[, 1],
        param_est = msm_param,
        ci_high = ci_msm_param[, 2],
        param_se = msm_se,
        p_value = pval_msm_param
      ) %>%
      as_tibble()
    }
  }
  # return results
  if (param_type == "shift") {
    results_this_iter <- bind_rows(results_shifted_means)
  } else if (param_type == "npmsm") {
    results_this_iter <- msm_out[2, -1]
  }
  return(results_this_iter)
}
