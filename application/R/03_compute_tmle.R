library(here)
library(data.table)
library(tidyverse)
library(hal9001)
library(haldensify)
library(txshift)
set.seed(78914)

# load nuisance parameter estimates for TMLE computation
cd4_est_data <- readRDS(file = here("data", "cd4_nuisance_delta_grid.rds"))
cd8_est_data <- readRDS(file = here("data", "cd8_nuisance_delta_grid.rds"))

# load updated Janes + Fong et al. data set and reduce for analysis
hvtn505_data <- read_csv(here("data", "primary505_for_sharing.csv")) %>%
  mutate(
    racecc = case_when(racecc == "White" ~ 0,
                       racecc == "Black" ~ 1,
                       racecc == "Hispanic/Other" ~ 2)
  ) %>%
  select(
    pub_id,                              # for censoring indicator
    age, racecc, BMI, bhvrisk,           # baseline covariates
    HIVwk28preunbl                       # indicator of failure
  )

# find missing values in baseline covariates
hvtn505_missing_baseline <-
  which(is.na(hvtn505_data$bhvrisk) | is.na(hvtn505_data$BMI))

janes_data <- read_csv(here("data", "janes_ipc_weights.csv")) %>%
  select(pub_id,                              # for censoring indicator
         age, racecc, BMI, bhvrisk,           # baseline covariates
         cd4_score_scaled, cd8_score_scaled,  # immune markers from Janes et al
         HIVwk28preunbl,                      # indicator of failure
         wt_hal                               # IPC weights based on HAL
        )

# get IPC weights for all study participants
hvtn505_ipcw <- read_csv(here("data", "hvtn505_ipc_weights.csv"))

# components of data needed for estimation
cens_ind <- (hvtn505_data$pub_id %in% janes_data$pub_id) %>%
  enframe(name = NULL) %>%
  slice(-hvtn505_missing_baseline) %>%
  unlist() %>%
  as.numeric()
w_baseline <- hvtn505_data %>%
  select(age, racecc, BMI, bhvrisk) %>%
  slice(-hvtn505_missing_baseline) %>%
  as.matrix()
y_surv_ind <- hvtn505_data %>%
  slice(-hvtn505_missing_baseline) %>%
  select(HIVwk28preunbl) %>%
  unlist() %>%
  as.numeric()
ipcw_estim <- list(pi_mech = hvtn505_ipcw$pi_mech,
                   ipc_weights = hvtn505_ipcw$ipc_weights[cens_ind == 1]
                  )
hvtn505_data <- hvtn505_data %>%
  slice(-hvtn505_missing_baseline)

################################################################################
## compute TML estimator via wrapper function
################################################################################

do_shift_tmle <- function(gn_estim_in, Qn_estim_in, A_in) {
  # set correct names for g and Q components for input
  gn_estim <- gn_estim_in %>%
    rename(downshift = gn_downshift, noshift = gn_natural,
           upshift = gn_upshift, upupshift = gn_upupshift) %>%
    as.data.table()

  Qn_estim <- Qn_estim_in %>%
    transmute(
      noshift = txshift:::scale_to_unit(Qn_natural),
      upshift = txshift:::scale_to_unit(Qn_shifted)
    ) %>%
    as.data.table()

  # compute TML estimator
  shift_tmle <- txshift(W = w_baseline,
                        A = A_in,
                        Y = y_surv_ind,
                        C = cens_ind,
                        V = c("W", "Y"),
                        max_iter = 1,
                        estimator = "tmle",
                        fluctuation = "standard",
                        ipcw_fit_args = list(fit_type = "fit_spec"),
                        g_fit_args = list(fit_type = "fit_spec"),
                        Q_fit_args = list(fit_type = "fit_spec"),
                        ipcw_fit_spec = ipcw_estim,
                        gn_fit_spec = gn_estim,
                        Qn_fit_spec = Qn_estim,
                        eif_reg_type = "hal")
  return(shift_tmle)
}

################################################################################
## compute TML estimator for shifted CD4 data
################################################################################

# create A vector with missing values for immune responses of censored group
cd4_interv <- rep(0, nrow(hvtn505_data))
cd4_interv[cens_ind == 1] <- janes_data$cd4_score_scaled

cd4_shifted <- lapply(seq_along(cd4_est_data$gn), function(x) {
  gn_estim_in <- cd4_est_data$gn[[x]]
  Qn_estim_in <- cd4_est_data$Qn[[x]]
  out <- do_shift_tmle(gn_estim_in, Qn_estim_in, A_in = cd4_interv)
  return(out)
})


################################################################################
## compute TML estimator for shifted CD8 data
################################################################################

# create A vector with missing values for immune responses of censored group
cd8_interv <- rep(0, nrow(hvtn505_data))
cd8_interv[cens_ind == 1] <- janes_data$cd8_score_scaled

cd8_shifted <- lapply(seq_along(cd8_est_data$gn), function(x) {
  gn_estim_in <- cd8_est_data$gn[[x]]
  Qn_estim_in <- cd8_est_data$Qn[[x]]
  out <- do_shift_tmle(gn_estim_in, Qn_estim_in, A_in = cd8_interv)
  return(out)
})


################################################################################
## save results of TMLE computation
################################################################################

# save results from TML estimation
tmle_out <- list(cd4_tmle_results = cd4_shifted,
                 cd8_tmle_results = cd8_shifted)
saveRDS(tmle_out, file = here("data", "tmle_results_hal_eifreg.rds"))
