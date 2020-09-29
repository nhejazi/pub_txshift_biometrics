library(here)
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggsci)
library(knitr)
options(scipen = 999)
pd <- position_dodge(0.5)

# confidence interval multiplier
ci_level <- 0.95
ci_mult <- abs(qnorm(p = (1 - ci_level) / 2))

# simulation type: set data subdirectory based on one of three types
sim_type <- "hvtn505_null"
source(here("R", "utils", "summary_helpers.R"))
source(here("R", "utils", "get_truth.R"))
estimators <- c("Efficient w/ HAL", "Efficient w/ GLM",
                "Inefficient w/ HAL", "Inefficient w/ GLM")

# compute truth
truths <- lapply(delta_shift, function(delta) {
  get_truth_hvtn505(delta = delta, type = unlist(str_split(sim_type, "_"))[2])
})
truths <- do.call(c, truths)
msm_truth <- coef(lm(truths ~ delta_shift))[2]

sink(here("tables", "null_truths_hvtn.tex"))
rbind(as.integer(delta_shift), truths) %>%
  kable(format = "latex", digits = 4, label = "table:hvtn_null_truths",
        caption = "True counterfactual means across posited shift values") %>%
  kable_styling()
sink()

################################################################################
# SUMMARIZE SHIFT PARAMETER ESTIMATES OVER GRID
################################################################################

# load data
## efficient TMLE with HAL for censoring fit
file_names_01_tmle <- str_subset(dir(here("data", sim_type)),
                                 "01-eff_tmle-pi_hal_shift_n")
file_names_01_tmle_recent <- get_recent_files(file_names_01_tmle, num_files = 1)
eff_tmle_pi_hal <- readRDS(here("data", sim_type, file_names_01_tmle_recent))

## efficient one-step with HAL for censoring fit
file_names_01_os <- str_subset(dir(here("data", sim_type)),
                               "01-eff_onestep-pi_hal_shift_n")
file_names_01_os_recent <- get_recent_files(file_names_01_os, num_files = 1)
eff_os_pi_hal <- readRDS(here("data", sim_type, file_names_01_os_recent))

# get samples sizes and re-order efficient simulations
n_obs <- as.numeric(str_remove(str_extract(file_names_01_tmle_recent,
                                           "n...."), "n"))

# helper function to create simulation summary statistics
summarise_simulation <- function(simulation_results) {
    simulation_summary <- simulation_results %>%
      bind_rows() %>%
      mutate(
        delta = round(delta, 1),
        bias = case_when(delta == delta_shift[1] ~ (truths[1] - param_est),
                         delta == delta_shift[2] ~ (truths[2] - param_est),
                         delta == delta_shift[3] ~ (truths[3] - param_est),
                         delta == delta_shift[4] ~ (truths[4] - param_est),
                         delta == delta_shift[5] ~ (truths[5] - param_est),
                         delta == delta_shift[6] ~ (truths[6] - param_est),
                         delta == delta_shift[7] ~ (truths[7] - param_est),
                         delta == delta_shift[8] ~ (truths[8] - param_est),
                         delta == delta_shift[9] ~ (truths[9] - param_est)),
        covers = case_when(delta == delta_shift[1] ~
                             truths[1] >= lwr_ci & truths[1] <= upr_ci,
                           delta == delta_shift[2] ~
                             truths[2] >= lwr_ci & truths[2] <= upr_ci,
                           delta == delta_shift[3] ~
                             truths[3] >= lwr_ci & truths[3] <= upr_ci,
                           delta == delta_shift[4] ~
                             truths[4] >= lwr_ci & truths[4] <= upr_ci,
                           delta == delta_shift[5] ~
                             truths[5] >= lwr_ci & truths[5] <= upr_ci,
                           delta == delta_shift[6] ~
                             truths[6] >= lwr_ci & truths[6] <= upr_ci,
                           delta == delta_shift[7] ~
                             truths[7] >= lwr_ci & truths[7] <= upr_ci,
                           delta == delta_shift[8] ~
                             truths[8] >= lwr_ci & truths[8] <= upr_ci,
                           delta == delta_shift[9] ~
                             truths[9] >= lwr_ci & truths[9] <= upr_ci)
      ) %>%
      group_by(delta) %>%
      summarise(
        est_avg = mean(param_est),
        var_mc = var(param_est),
        var_avg = mean(param_var),
        coverage = mean(covers),
        cov_var_mc = coverage * (1 - coverage),
        bias_avg = mean(bias),
        bias_var_mc = var(bias),
        mse_est = mean(bias^2 + var_mc),
        mse_var_mc = var(bias^2 + var_mc),
        n_sim = n()
      ) %>%
      ungroup() %>%
      mutate(
        bias_avg_ci_lwr = bias_avg - sqrt(bias_var_mc / n_sim) * ci_mult,
        bias_avg_ci_upr = bias_avg + sqrt(bias_var_mc / n_sim) * ci_mult,
        mse_est_ci_lwr = mse_est - sqrt(mse_var_mc / n_sim) * ci_mult,
        mse_est_ci_upr = mse_est + sqrt(mse_var_mc / n_sim) * ci_mult,
        cov_est_ci_lwr = coverage - sqrt(cov_var_mc / n_sim) * ci_mult,
        cov_est_ci_upr = coverage + sqrt(cov_var_mc / n_sim) * ci_mult
      )
  return(simulation_summary)
}

# create simulation summaries
eff_tmle_pi_hal_summary <- summarise_simulation(eff_tmle_pi_hal)
eff_os_pi_hal_summary <- summarise_simulation(eff_os_pi_hal)

# put all results together for easier plotting
sim_est_results <- list(eff_tmle_pi_hal = eff_tmle_pi_hal_summary,
                        eff_os_pi_hal = eff_os_pi_hal_summary) %>%
  bind_rows(.id = "est_type") %>%
  mutate(
    estimator = case_when(est_type == "eff_tmle_pi_hal" ~ "TMLE",
                          est_type == "eff_os_pi_hal" ~ "One-step")
  )

# plot of absolute bias
p_bias <- sim_est_results %>%
  ggplot(aes(x = as.factor(delta), y = bias_avg, group = estimator,
             shape = estimator, fill = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = bias_avg_ci_lwr, ymax = bias_avg_ci_upr),
                position = pd, width = 0.25, linetype = "dotted") +
  geom_line(linetype = "dashed", position = pd) +
  geom_point(size = 6, alpha = 0.75, position = pd) +
  ylim(-0.0025, 0.0005) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background =
          element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
        legend.title = element_blank(),
        text = element_text(size = 25),
        axis.text.x = element_text(colour = "black", size = 22, hjust = 1,
                                   angle = 30),
        axis.text.y = element_text(colour = "black", size = 22)) +
  labs(x = TeX("Shift value $\\delta$"),
       y = TeX("$\\psi - \\hat{\\psi}$"),
       title = "Average bias"
      ) +
  scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
  scale_fill_nejm()
ggsave(filename = here("graphs", "manuscript",
                       paste(sim_type, "bias_avg.pdf", sep = "_")),
       plot = p_bias, width = 19, height = 13)

# plot of MSE
p_mse <- sim_est_results %>%
  ggplot(aes(x = as.factor(delta), y = mse_est, group = estimator,
             shape = estimator, fill = estimator)) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linetype = "dashed", position = pd) +
  geom_errorbar(aes(ymin = mse_est_ci_lwr, ymax = mse_est_ci_upr),
                position = pd, width = 0.25, linetype = "dotted") +
  geom_point(size = 6, alpha = 0.75, position = pd) +
  ylim(0.00007, 0.000085) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background =
          element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
        legend.title = element_blank(),
        text = element_text(size = 25),
        axis.text.x = element_text(colour = "black", size = 22, hjust = 1,
                                   angle = 30),
        axis.text.y = element_text(colour = "black", size = 22)) +
  labs(x = TeX("Shift value $\\delta$"),
       y = TeX("MSE = $(\\psi - \\hat{\\psi})^2 + \\hat{\\sigma}^2$"),
       title = "Mean squared error"
      ) +
  scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
  scale_fill_nejm()
ggsave(filename = here("graphs", "manuscript",
                       paste(sim_type, "mse.pdf", sep = "_")),
       plot = p_mse, width = 19, height = 13)

# plot of coverage
p_cov <- sim_est_results %>%
  ggplot(aes(x = as.factor(delta), y = coverage, group = estimator,
             shape = estimator, fill = estimator)) +
  geom_hline(yintercept = ci_level, linetype = "dashed") +
  geom_line(linetype = "dashed", position = pd) +
  geom_errorbar(aes(ymin = cov_est_ci_lwr, ymax = cov_est_ci_upr),
                position = pd, width = 0.25, linetype = "dotted") +
  geom_point(size = 6, alpha = 0.75, position = pd) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.background =
          element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
        legend.title = element_blank(),
        text = element_text(size = 25),
        axis.text.x = element_text(colour = "black", size = 22, hjust = 1,
                                   angle = 30),
        axis.text.y = element_text(colour = "black", size = 22)) +
  labs(x = TeX("Shift value $\\delta$"),
       y = "Empirical coverage",
       title = "Coverage of 95% CIs"
      ) +
  scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
  scale_fill_nejm()
ggsave(filename = here("graphs", "manuscript",
                       paste(sim_type, "cov.pdf", sep = "_")),
       plot = p_cov, width = 19, height = 13)

# panel plot of bias and MSE
p_bias <- p_bias + theme(legend.position = "none")
p_cov <- p_cov + theme(legend.position = "none")
p_panel <- p_bias + p_mse + p_cov
ggsave(filename = here("graphs", "manuscript",
                       paste(sim_type, "panel.pdf", sep = "_")),
       plot = p_panel, width = 19, height = 7)
