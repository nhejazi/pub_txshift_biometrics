# programmatic housekeeping
library(here)
library(tidyverse)
library(hal9001)

# load Janes + Fong et al. data set as well as full HVTN 505 data
janes_data <- read_csv(here("data", "hvtn505_for_analysis.csv"))
primary_data <- read_csv(here("data", "primary505_for_sharing.csv"))
delta_ind <- primary_data %>%
  dplyr::filter(!is.na(bhvrisk) & !is.na(BMI)) %>%
  transmute(
    delta = as.numeric(pub_id %in% janes_data$pub_id)
  ) %>%
  unlist() %>%
  as.numeric()

# determine covariates informative of sampling mechanism based on overlap
# between Janes et al. and primary HVTN 505 data
v_names <- janes_data[, colnames(janes_data) %in% colnames(primary_data)] %>%
    select(-pub_id) %>%
    colnames()

# Janes et al. use treatment group, BMI, and race/ethnicity
v_covars <- primary_data %>%
    select(v_names) %>%
    mutate(
      racecc = case_when(racecc == "White" ~ 0,
                         racecc == "Black" ~ 1,
                         !(racecc %in% c("White", "Black")) ~ 2
                        ),
      hg_strata = as.numeric(as.factor(hg_strata))
    ) %>%
    select(HIVpreunbl, BMI, racecc, age, bhvrisk) %>%
    dplyr::filter(!is.na(bhvrisk) & !is.na(BMI))

# compute IPC weights using HAL
ipc_hal_fit <- fit_hal(X = v_covars, Y = delta_ind,
                       fit_type = "glmnet",
                       family = "binomial",
                       standardize = FALSE,
                       lambda.min.ratio = 0.00001)
ipc_hal_pred <- predict(ipc_hal_fit, new_data = v_covars)
ipc_weights <- delta_ind / ipc_hal_pred
ipcw_out <- as_tibble(list(pi_mech = ipc_hal_pred, ipc_weights = ipc_weights))

# add HAL IPC weights to Janes et al. data
janes_data  <- janes_data %>%
  mutate(
    wt_hal = ipc_weights[delta_ind == 1],
    racecc = case_when(racecc == "White" ~ 0,
                       racecc == "Black" ~ 1,
                       !(racecc %in% c("White", "Black")) ~ 2
                      ),
    hg_strata = as.numeric(as.factor(hg_strata))
  )

# output intermediary data
write_csv(x = janes_data, path = here("data", "janes_ipc_weights.csv"))
write_csv(x = ipcw_out, path = here("data", "hvtn505_ipc_weights.csv"))

