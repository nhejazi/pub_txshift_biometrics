msm_shift <- function(tmle_fit_estimates,
                      delta_grid,
                      msm_form = list(type = "linear", knot = NA),
                      level = 0.95,
                      weights = NULL) {

  # make sure more than one parameter has been estimated for trend
  assertthat::assert_that(length(tmle_fit_estimates) > 1)

  # matrix of EIF(O_i) values and estimates across each parameter estimated
  eif_mat <- sapply(tmle_fit_estimates, `[[`, "eif")
  psi_vec <- sapply(tmle_fit_estimates, `[[`, "psi")

  # set weights to be the inverse of the variance of each TML estimate
  if (is.null(weights)) {
    weights <- as.numeric(1 / diag(stats::cov(eif_mat)))
  }

  # multiplier for CI construction
  ci_mult <- (c(1, -1) * stats::qnorm((1 - level) / 2))

  # set right-hand side of MSM formula
  if (msm_form[["type"]] == "piecewise") {
    msm_rhs <- paste0("delta + I(pmax(delta - ", msm_form[["knot"]], ", 0))")
  } else {
    msm_rhs <- "delta"
  }

  # create design matrix for MSM
  x_mat <- stats::model.matrix(
    stats::as.formula(paste0("psi_vec ~ ", msm_rhs)),
    data = data.frame(psi_vec = psi_vec, delta = delta_grid)
  )

  # compute the MSM parameters
  omega <- diag(weights)
  s_mat <- solve(t(x_mat) %*% omega %*% x_mat) %*% t(x_mat) %*% omega
  msm_param <- as.vector(s_mat %*% psi_vec)

  # compute inference for MSM based on individual EIF(O_i) for each parameter
  msm_eif <- tcrossprod(eif_mat, s_mat)
  msm_var <- diag(stats::cov(msm_eif))
  msm_se <- sqrt(msm_var / nrow(msm_eif))

  # build confidence intervals and hypothesis tests for EIF(msm)
  ci_msm_param <- msm_se %*% t(ci_mult) + msm_param
  pval_msm_param <- 2 * stats::pnorm(-abs(msm_param / msm_se))

  # tables for output
  txshift_out <- list(
    ci_low = psi_vec + ci_mult[1] * sqrt(diag(stats::cov(eif_mat))),
    psi = psi_vec,
    ci_high = psi_vec + ci_mult[2] * sqrt(diag(stats::cov(eif_mat)))
  ) %>%
  tibble::as_tibble()

  msm_out <- list(
    param = names(msm_se),
    ci_low = ci_msm_param[, 1],
    param_est = msm_param,
    ci_high = ci_msm_param[, 2],
    param_se = msm_se,
    p_value = pval_msm_param
  ) %>%
  tibble::as_tibble()

  out <- list(marginal = txshift_out, msm = msm_out)
  return(out)
}
