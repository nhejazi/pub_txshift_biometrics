library(here)
library(tidyverse)
library(latex2exp)
library(patchwork)
library(ggsci)
library(knitr)
options(scipen = 999)

# helper functions for ggplot2
pd <- position_dodge(0.5)
label_maker <- c(
  "-0.5" = "Shift down (-1/2 std. unit)",
  "0" = "No shift (0 std. units)",
  "0.5" = "Shift up (+1/2 std. unit)",
  "ipcw_eif" = "Augmented EIF (efficient with ML)",
  "ipcw_loss" = "Augmented Loss (efficient with GLMs)",
  "gn_hal" = "HAL for sampling mechanism",
  "gn_glm" = "GLM for sampling mechanism"
)

# helper function to get 5 most recent files for efficient estimator simulations
get_recent_files <- function(file_names, num_files = 10) {
  file_dates <- as.Date(str_extract(file_names, "2019......"))
  recent_files <- sort(file_dates, decreasing = TRUE)[seq_len(num_files)]
  recent_files_regex <- paste(recent_files, collapse = '|')
  which_recent <- str_detect(file_names, recent_files_regex)
  return(file_names[which_recent])
}

# load data
## efficient TMLE with HAL for censoring fit
file_names_01_tmle <- str_subset(dir(here("data", "simple")),
                                 "01-eff_tmle-pi_hal_n")
file_names_01_tmle_recent <- get_recent_files(file_names_01_tmle, num_files = 5)
eff_tmle_pi_hal <- lapply(file_names_01_tmle_recent, function(file) {
  readRDS(here("data", "simple", file))
})
names(eff_tmle_pi_hal) <- file_names_01_tmle_recent %>%
  str_extract(regex("n....")) %>%
  str_remove("_") %>%
  str_replace("n", "n_")

## efficient TMLE with GLM for censoring fit
file_names_02_tmle <- str_subset(dir(here("data", "simple")),
                                 "02-eff_tmle-pi_glm_n")
file_names_02_tmle_recent <- get_recent_files(file_names_02_tmle, num_files = 5)
eff_tmle_pi_glm <- lapply(file_names_02_tmle_recent, function(file) {
  readRDS(here("data", "simple", file))
})
names(eff_tmle_pi_glm) <- file_names_02_tmle_recent %>%
  str_extract(regex("n....")) %>%
  str_remove("_") %>%
  str_replace("n", "n_")

## inefficient TMLE with HAL for censoring fit
file_names_03_tmle <- str_subset(dir(here("data", "simple")),
                                 "03-ineff_tmle-pi_hal")
file_names_03_tmle_recent <- get_recent_files(file_names_03_tmle, num_files = 1)
ineff_tmle_pi_hal <- readRDS(here("data", "simple", file_names_03_tmle_recent))

## inefficient TMLE with GLM for censoring fit
file_names_04_tmle <- str_subset(dir(here("data", "simple")),
                                 "04-ineff_tmle-pi_glm")
file_names_04_tmle_recent <- get_recent_files(file_names_04_tmle, num_files = 1)
ineff_tmle_pi_glm <- readRDS(here("data", "simple", file_names_04_tmle_recent))

## inefficient TMLE ignoring 2nd-stage sampling
file_names_05_tmle <- str_subset(dir(here("data", "simple")),
                                 "05-ignore_sampling-tmle")
file_names_05_tmle_recent <- get_recent_files(file_names_05_tmle, num_files = 1)
ignore_sampling_tmle <- readRDS(here("data", "simple",
                                     file_names_05_tmle_recent))

## efficient one-step with HAL for censoring fit
file_names_01_os <- str_subset(dir(here("data", "simple")),
                               "01-eff_onestep-pi_hal_n")
file_names_01_os_recent <- get_recent_files(file_names_01_os, num_files = 5)
eff_os_pi_hal <- lapply(file_names_01_os_recent, function(file) {
  readRDS(here("data", "simple", file))
})
names(eff_os_pi_hal) <- file_names_01_os_recent %>%
  str_extract(regex("_n....")) %>%
  str_remove_all("_") %>%
  str_replace("n", "n_")

## efficient one-step with GLM for censoring fit
file_names_02_os <- str_subset(dir(here("data", "simple")),
                               "02-eff_onestep-pi_glm_n")
file_names_02_os_recent <- get_recent_files(file_names_02_os, num_files = 5)
eff_os_pi_glm <- lapply(file_names_02_os_recent, function(file) {
  readRDS(here("data", "simple", file))
})
names(eff_os_pi_glm) <- file_names_02_os_recent %>%
  str_extract(regex("_n....")) %>%
  str_remove_all("_") %>%
  str_replace("n", "n_")

## inefficient one-step with HAL for censoring fit
file_names_03_os <- str_subset(dir(here("data", "simple")),
                               "03-ineff_onestep-pi_hal")
file_names_03_os_recent <- get_recent_files(file_names_03_os, num_files = 1)
ineff_os_pi_hal <- readRDS(here("data", "simple", file_names_03_os_recent))

## inefficient one-step with GLM for censoring fit
file_names_04_os <- str_subset(dir(here("data", "simple")),
                               "04-ineff_onestep-pi_glm")
file_names_04_os_recent <- get_recent_files(file_names_04_os, num_files = 1)
ineff_os_pi_glm <- readRDS(here("data", "simple", file_names_04_os_recent))

## inefficient one-step ignoring 2nd-stage sampling
file_names_05_os <- str_subset(dir(here("data", "simple")),
                               "05-ignore_sampling-onestep")
file_names_05_os_recent <- get_recent_files(file_names_05_os, num_files = 1)
ignore_sampling_os <- readRDS(here("data", "simple", file_names_05_os_recent))

# get samples sizes and re-order efficient simulations
estimators <- c("Efficient w/ HAL", "Efficient w/ GLM",
                "Inefficient w/ HAL", "Inefficient w/ GLM")
n_obs <- as.numeric(str_remove(names(eff_tmle_pi_hal), "n_"))
n_obs_idx <- sort.int(n_obs, index.return = TRUE)$ix
eff_tmle_pi_hal <- eff_tmle_pi_hal[n_obs_idx]
eff_tmle_pi_glm <- eff_tmle_pi_glm[n_obs_idx]
eff_os_pi_hal <- eff_os_pi_hal[n_obs_idx]
eff_os_pi_glm <- eff_os_pi_glm[n_obs_idx]

# get truth for simulation
source(here("R", "get_truth.R"))
truth_downshift <- get_truth(delta = -0.5)
truth_noshift <- get_truth(delta = 0)
truth_upshift <- get_truth(delta = 0.5)


# helper function to create simulation summary statistics
summarise_simulation <- function(simulation_results) {
  simulation_summary <- bind_rows(simulation_results, .id = "n_samp") %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n_")),
      bias = case_when(delta == -0.5 ~ (param_est - truth_downshift),
                       delta == 0 ~ (param_est - truth_noshift),
                       delta == 0.5 ~ (param_est - truth_upshift)),
      covers = case_when(delta == -0.5 ~ data.table::between(truth_downshift,
                                                             lwr_ci, upr_ci),
                         delta == 0 ~ data.table::between(truth_noshift,
                                                          lwr_ci, upr_ci),
                         delta == 0.5 ~ data.table::between(truth_upshift,
                                                            lwr_ci, upr_ci))
    ) %>%
    group_by(delta, n_samp) %>%
    summarise(
      est_avg = mean(param_est),
      var_mc = var(param_est),
      var_avg = mean(param_var),
      coverage = mean(covers),
      bias_abs = abs(mean(bias)),
      bias_sqrtn = bias_abs * sqrt(unique(n_samp)),
      nmse = (mean(bias)^2 + var_mc) * unique(n_samp),
      n_sim = n()
    ) %>%
    ungroup()
  return(simulation_summary)
}

# create simulation summaries
eff_tmle_pi_hal_summary <- summarise_simulation(eff_tmle_pi_hal)
eff_tmle_pi_glm_summary <- summarise_simulation(eff_tmle_pi_glm)
ineff_tmle_pi_hal_summary <- summarise_simulation(ineff_tmle_pi_hal)
ineff_tmle_pi_glm_summary <- summarise_simulation(ineff_tmle_pi_glm)
naive_tmle_ignore_summary <- summarise_simulation(ignore_sampling_tmle)
eff_os_pi_hal_summary <- summarise_simulation(eff_os_pi_hal)
eff_os_pi_glm_summary <- summarise_simulation(eff_os_pi_glm)
ineff_os_pi_hal_summary <- summarise_simulation(ineff_os_pi_hal)
ineff_os_pi_glm_summary <- summarise_simulation(ineff_os_pi_glm)
naive_os_ignore_summary <- summarise_simulation(ignore_sampling_os)

# put all results together for easier plotting
est_results <- list(eff_tmle_pi_hal = eff_tmle_pi_hal_summary,
                    eff_tmle_pi_glm = eff_tmle_pi_glm_summary,
                    ineff_tmle_pi_hal = ineff_tmle_pi_hal_summary,
                    ineff_tmle_pi_glm = ineff_tmle_pi_glm_summary,
                    naive_tmle = naive_tmle_ignore_summary,
                    eff_os_pi_hal = eff_os_pi_hal_summary,
                    eff_os_pi_glm = eff_os_pi_glm_summary,
                    ineff_os_pi_hal = ineff_os_pi_hal_summary,
                    ineff_os_pi_glm = ineff_os_pi_glm_summary,
                    naive_os = naive_os_ignore_summary
                   )
sim_est_results <- bind_rows(est_results, .id = "est_type") %>%
  mutate(
    estimator = case_when(est_type == "eff_tmle_pi_hal" ~ "Augmented TMLE",
                          est_type == "eff_tmle_pi_glm" ~ "Augmented TMLE",
                          est_type == "ineff_tmle_pi_hal" ~ "Reweighted TMLE",
                          est_type == "ineff_tmle_pi_glm" ~ "Reweighted TMLE",
                          est_type == "naive_tmle" ~ "Naive TMLE",
                          est_type == "eff_os_pi_hal" ~ "Augmented one-step",
                          est_type == "eff_os_pi_glm" ~ "Augmented one-step",
                          est_type == "ineff_os_pi_hal" ~ "Reweighted one-step",
                          est_type == "ineff_os_pi_glm" ~ "Reweighted one-step",
                          est_type == "naive_os" ~ "Naive one-step"
                         ),
    correction = case_when(est_type == "eff_tmle_pi_hal" ~ "ipcw_eif",
                           est_type == "eff_tmle_pi_glm" ~ "ipcw_eif",
                           est_type == "ineff_tmle_pi_hal" ~ "ipcw_loss",
                           est_type == "ineff_tmle_pi_glm" ~ "ipcw_loss",
                           est_type == "naive_tmle" ~ "ipcw_loss",
                           est_type == "eff_os_pi_hal" ~ "ipcw_eif",
                           est_type == "eff_os_pi_glm" ~ "ipcw_eif",
                           est_type == "ineff_os_pi_hal" ~ "ipcw_loss",
                           est_type == "ineff_os_pi_glm" ~ "ipcw_loss",
                           est_type == "naive_os" ~ "ipcw_loss"
                          ),
    gn_type = case_when(est_type == "eff_tmle_pi_hal" ~ "gn_hal",
                        est_type == "eff_tmle_pi_glm" ~ "gn_glm",
                        est_type == "ineff_tmle_pi_hal" ~ "gn_hal",
                        est_type == "ineff_tmle_pi_glm" ~ "gn_glm",
                        est_type == "naive_tmle" ~ "gn_glm",
                        est_type == "eff_os_pi_hal" ~ "gn_hal",
                        est_type == "eff_os_pi_glm" ~ "gn_glm",
                        est_type == "ineff_os_pi_hal" ~ "gn_hal",
                        est_type == "ineff_os_pi_glm" ~ "gn_glm",
                        est_type == "naive_os" ~ "gn_glm"
                       )
  )


for (delta_shift in unique(sim_est_results$delta)) {
  # plot of root-n scaled bias
  p_nbias <- sim_est_results %>%
    dplyr::filter(delta == delta_shift) %>%
    ggplot(aes(x = as.factor(n_samp), y = bias_sqrtn, group = estimator,
               shape = estimator, fill = estimator)) +
    geom_point(size = 8, position = pd) +
    geom_line(linetype = "dotted", position = pd) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none",
          legend.background =
            element_rect(fill = "gray90", size = 0.25, linetype = "dotted"),
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(colour = "black", size = 22),
          axis.text.y = element_text(colour = "black", size = 22)) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "Sample size",
         y = TeX("$\\sqrt{n} \\times$ |$\\psi - \\hat{\\psi}$|"),
         title = TeX(paste("Scaled bias,", "$\\delta$ =", delta_shift))
        ) +
    scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
    scale_fill_nejm() +
    facet_grid(vars(gn_type), labeller = as_labeller(label_maker))
  ggsave(filename = here("graphs", "manuscript",
                         paste0("simple_effect_nbias_delta_",
                                delta_shift, ".pdf")),
         plot = p_nbias, width = 19, height = 13)

  # plot of n-scaled MSE
  p_nmse <- sim_est_results %>%
    dplyr::filter(delta == delta_shift) %>%
    ggplot(aes(x = as.factor(n_samp), y = nmse, group = estimator,
               shape = estimator, fill = estimator)) +
    geom_point(size = 8, position = pd) +
    geom_line(linetype = "dotted", position = pd) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.background =
            element_rect(fill = "gray90", size = .25, linetype = "dotted"),
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(colour = "black", size = 22),
          axis.text.y = element_text(colour = "black", size = 22)) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "Sample size",
         y = TeX(paste("$n \\cdot MSE = n \\times",
                       "((\\psi - \\hat{\\psi})^2 + \\hat{\\sigma}^2)$")),
         title = TeX(paste("Scaled MSE,",
                           "$\\delta$ =", delta_shift))
        ) +
    scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
    scale_fill_nejm() +
    facet_grid(vars(gn_type), labeller = as_labeller(label_maker),
               scales = "free_y")
  ggsave(filename = here("graphs", "manuscript",
                         paste0("simple_effect_nmse_delta_",
                                delta_shift, ".pdf")),
         plot = p_nmse, width = 19, height = 13)

  # plot of confidence interval coverage
  p_cov <- sim_est_results %>%
    dplyr::filter(delta == delta_shift) %>%
    ggplot(aes(x = as.factor(n_samp), y = coverage, group = estimator,
               shape = estimator, fill = estimator)) +
    geom_point(size = 8, position = pd) +
    geom_line(linetype = "dotted", position = pd) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "none",
          legend.background =
            element_rect(fill = "gray90", size = .25, linetype = "dotted"),
          legend.title = element_blank(),
          text = element_text(size = 25),
          axis.text.x = element_text(colour = "black", size = 22),
          axis.text.y = element_text(colour = "black", size = 22)) +
    guides(color = guide_legend(title = NULL)) +
    labs(x = "Sample size",
         y = paste("Coverage across", unique(sim_est_results$n_sim),
                   "simulations"),
         title = TeX(paste("CI coverage,",
                           "$\\delta$ =", delta_shift))
        ) +
    scale_shape_manual(values = c(22, 25, 12, 13, 23, 24)) +
    scale_fill_nejm() +
    facet_grid(vars(gn_type), labeller = as_labeller(label_maker))
  ggsave(filename = here("graphs", "manuscript",
                         paste0("simple_effect_coverage_delta_",
                                delta_shift, ".pdf")),
         plot = p_cov, width = 19, height = 13)

  # now, panel plot for specified shift
  p_panel <- p_nbias + p_nmse + p_cov
  ggsave(filename = here("graphs", "manuscript",
                         paste0("simple_effect_panel_delta_",
                                delta_shift, ".pdf")),
         plot = p_panel, width = 19, height = 13)
}
