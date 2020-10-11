library(here)
library(data.table)
library(tidyverse)
library(patchwork)
library(latex2exp)
library(hal9001)
library(haldensify)
library(lspline)
library(txshift)
set.seed(78914)
source(here("R", "trend_test.R"))

# for use with scaled CD4/CD8 score
delta_shift_grid <- seq(-2, 2, by = 0.5)

# load TMLE results for summarization
tmle_shift_results <- readRDS(file = here("data",
                                          "tmle_results_hal_eifreg.rds"))
cd4_results <- tmle_shift_results$cd4_tmle_results
cd8_results <- tmle_shift_results$cd8_tmle_results


# fit linear and spline MSMs for grid of shifts
cd4_msm_linear <- msm_shift(tmle_fit_estimates = cd4_results,
                            delta_grid = delta_shift_grid,
                            weights = rep(1, length(cd4_results)))
cd4_msm_spline <- msm_shift(tmle_fit_estimates = cd4_results,
                            delta_grid = delta_shift_grid,
                            msm_form = list(type = "piecewise", knot = 0),
                            weights = rep(1, length(cd4_results)))
cd8_msm_linear <- msm_shift(tmle_fit_estimates = cd8_results,
                            delta_grid = delta_shift_grid,
                            weights = rep(1, length(cd4_results)))
cd8_msm_spline <- msm_shift(tmle_fit_estimates = cd8_results,
                            delta_grid = delta_shift_grid,
                            msm_form = list(type = "piecewise", knot = 0),
                            weights = rep(1, length(cd4_results)))

# hacky stuff to get ggplot2 to also include the spline model
cd4_msm_data <- copy(cd4_msm_spline$marginal) %>%
  as.data.table() %>%
  set(j = "delta", value = delta_shift_grid) %>%
  setnames(c("psi", "delta"), c("y", "x")) %>%
  select(c("y", "x"))
cd4_msm_mod <- lm(y ~ lspline(x, 0, marginal = TRUE), data = cd4_msm_data)
cd8_msm_data <- copy(cd8_msm_spline$marginal) %>%
  as.data.table() %>%
  set(j = "delta", value = delta_shift_grid) %>%
  setnames(c("psi", "delta"), c("y", "x")) %>%
  select(c("y", "x"))
cd8_msm_mod <- lm(y ~ lspline(x, 0, marginal = TRUE), data = cd8_msm_data)

# plot results of TMLEs and MSM for CD4+
p_cd4_tmle <- sapply(cd4_results, confint) %>%
  t() %>%
  as_tibble() %>%
  add_column(delta = delta_shift_grid) %>%
  ggplot(aes(x = delta, y = est)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = lwr_ci, ymax = upr_ci),
                position = "dodge", linetype = "dotted", width = 0.2) +
  geom_segment(aes(x = min(delta_shift_grid), xend = max(delta_shift_grid),
                   y = cd4_msm_linear$msm$param_est[1] +
                     min(delta_shift_grid) * cd4_msm_linear$msm$param_est[2],
                   yend = cd4_msm_linear$msm$param_est[1] +
                    max(delta_shift_grid) * cd4_msm_linear$msm$param_est[2]),
               size = 0.5, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", formula = formula(cd4_msm_mod), se = FALSE,
              color = "black", size = 0.5, linetype = "dashed") +
  labs(
    x = "Posited change in standardized CD4+ polyfunctionality (sd units)",
    y = "Risk of HIV-1 infection in vaccine recipients",
    title = paste("TML estimates of mean counterfactual HIV-1 infection risk",
                  "under shifted CD4+ polyfunctionality"),
    subtitle = TeX(paste("with pointwise confidence intervals, working linear",
                         "MSM", paste0("($\\hat{\\beta}_{TMLE}$ =",
                         round(cd4_msm_linear$msm$param_est[2], 4), ")"),
                         ", and post-hoc spline MSM"))
  ) +
  theme_bw() +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(colour = "black", size = 32),
        axis.text.y = element_text(colour = "black", size = 32)) +
  scale_y_continuous(breaks = round(seq(0, 0.12, by = 0.02), 2),
                     limits = c(0, 0.11))
ggsave(filename = here::here("graphs", paste0("cd4_msm_tmle_summary.pdf")),
       plot = p_cd4_tmle, height = 13, width = 19)


# plot results of TMLEs and MSM for CD8+
p_cd8_tmle <- sapply(cd8_results, confint) %>%
  t() %>%
  as_tibble() %>%
  add_column(delta = delta_shift_grid) %>%
  ggplot(aes(x = delta, y = est)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = lwr_ci, ymax = upr_ci),
                position = "dodge", linetype = "dotted", width = 0.2) +
  geom_segment(aes(x = min(delta_shift_grid), xend = max(delta_shift_grid),
                   y = cd8_msm_linear$msm$param_est[1] +
                     min(delta_shift_grid) * cd8_msm_linear$msm$param_est[2],
                   yend = cd8_msm_linear$msm$param_est[1] +
                    max(delta_shift_grid) * cd8_msm_linear$msm$param_est[2]),
               size = 0.5, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", formula = formula(cd8_msm_mod), se = FALSE,
              color = "black", size = 0.5, linetype = "dashed") +
  labs(
    x = "Posited change in standardized CD8+ polyfunctionality (sd units)",
    y = "Risk of HIV-1 infection in vaccine recipients",
    title = paste("TML estimates of mean counterfactual HIV-1 infection risk",
                  "under shifted CD8+ polyfunctionality"),
    subtitle = TeX(paste("with pointwise confidence intervals, working linear",
                         "MSM", paste0("($\\hat{\\beta}_{TMLE}$ =",
                         round(cd8_msm_linear$msm$param_est[2], 4), ")"),
                         ", and post-hoc spline MSM"))
  ) +
  theme_bw() +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(colour = 'black', size = 32),
        axis.text.y = element_text(colour = 'black', size = 32)) +
  scale_y_continuous(breaks = round(seq(0, 0.12, by = 0.02), 2),
                     limits = c(0, 0.11))
ggsave(filename = here::here("graphs", paste0("cd8_msm_tmle_summary.pdf")),
       plot = p_cd8_tmle, height = 13, width = 19)

# combine plots
p_cdboth_tmle <- p_cd4_tmle / p_cd8_tmle
ggsave(filename = here::here("graphs", paste0("msm_tmle_summary_both.pdf")),
       plot = p_cdboth_tmle, height = 18, width = 26)
