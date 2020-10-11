library(here)
library(future)
library(tidyverse)
library(data.table)
library(hal9001)
library(haldensify)
library(SuperLearner)
library(sl3)
set.seed(78914)
plan(multiprocess)
sl3_debug_mode()

# propensity score and shifts with HAL-based conditional density estimation
delta_shift_grid <- seq(-2, 2, by = 0.5)  # for use with scaled CD4/CD8 score


# load updated Janes + Fong et al. data set and reduce for analysis
janes_data <- read_csv(here("data", "janes_ipc_weights.csv")) %>%
  select(age, racecc, BMI, bhvrisk,           # baseline covariates
         cd4_score_scaled, cd8_score_scaled,  # immune markers from Janes et al
         HIVwk28preunbl,                      # failure indicator
         wt_hal                               # IPC weights based on HAL
        )
w_baseline <- janes_data %>%
  select(age, racecc, BMI, bhvrisk) %>%
  as.matrix()


# first, for the CD4 data
gn_cd4_hal_model <- haldensify(W = as.matrix(w_baseline),
                               A = janes_data$cd4_score_scaled,
                               wts = janes_data$wt_hal)

gn_pred_hal_cd4 <- lapply(delta_shift_grid, function(delta_shift) {
  # predict for natural density
  gn_cd4_hal_natural <- predict(gn_cd4_hal_model,
                                new_W = as.matrix(w_baseline),
                                new_A = janes_data$cd4_score_scaled)

  # predict for upshifted density
  gn_cd4_hal_upshift <- predict(gn_cd4_hal_model,
                                new_W = as.matrix(w_baseline),
                                new_A = janes_data$cd4_score_scaled +
                                  delta_shift)

  # predict for 2x upshifted density (for bounds)
  gn_cd4_hal_upupshift <- predict(gn_cd4_hal_model,
                                  new_W = as.matrix(w_baseline),
                                  new_A = janes_data$cd4_score_scaled +
                                    2 * delta_shift)

  # predict for downshifted density
  gn_cd4_hal_downshift <- predict(gn_cd4_hal_model,
                                  new_W = as.matrix(w_baseline),
                                  new_A = janes_data$cd4_score_scaled -
                                    delta_shift)

  # reformat into single data frame for output
  gn_cd4_hal <- list(
      gn_downshift = gn_cd4_hal_downshift,
      gn_natural = gn_cd4_hal_natural,
      gn_upshift = gn_cd4_hal_upshift,
      gn_upupshift = gn_cd4_hal_upupshift
    ) %>%
    as_tibble()

  # output
  return(gn_cd4_hal)
})
names(gn_pred_hal_cd4) <- as.character(delta_shift_grid)


################################################################################
# now, for the CD8 data
################################################################################
gn_cd8_hal_model <- haldensify(W = as.matrix(w_baseline),
                               A = janes_data$cd8_score_scaled,
                               wts = janes_data$wt_hal)

gn_pred_hal_cd8 <- lapply(delta_shift_grid, function(delta_shift) {
  # predict for natural density
  gn_cd8_hal_natural <- predict(gn_cd8_hal_model,
                                new_W = as.matrix(w_baseline),
                                new_A = janes_data$cd8_score_scaled)

  # predict for upshifted density
  gn_cd8_hal_upshift <- predict(gn_cd8_hal_model,
                                new_W = as.matrix(w_baseline),
                                new_A = janes_data$cd8_score_scaled +
                                  delta_shift)

  # predict for 2x upshifted density (for bounds)
  gn_cd8_hal_upupshift <- predict(gn_cd8_hal_model,
                                  new_W = as.matrix(w_baseline),
                                  new_A = janes_data$cd8_score_scaled +
                                    2 * delta_shift)

  # predict for downshifted density
  gn_cd8_hal_downshift <- predict(gn_cd8_hal_model,
                                  new_W = as.matrix(w_baseline),
                                  new_A = janes_data$cd8_score_scaled -
                                          delta_shift)

  # reformat into single data frame for output
  gn_cd8_hal <- list(
      gn_downshift = gn_cd8_hal_downshift,
      gn_natural = gn_cd8_hal_natural,
      gn_upshift = gn_cd8_hal_upshift,
      gn_upupshift = gn_cd8_hal_upupshift
    ) %>%
    as_tibble()

  # output
  return(gn_cd8_hal)
})
names(gn_pred_hal_cd8) <- as.character(delta_shift_grid)


################################################################################
# compute outcome regression with SL including HAL
################################################################################
mean_lrnr <- Lrnr_mean$new()
glm_fast_lrnr <- Lrnr_glm_fast$new(family = binomial())
ridge_lrnr <- Lrnr_glmnet$new(alpha = 0, nfolds = 3)
lasso_lrnr <- Lrnr_glmnet$new(alpha = 1, nfolds = 3)
enet_lrnr_reg25 <- Lrnr_glmnet$new(alpha = 0.25, nfolds = 3)
enet_lrnr_reg50 <- Lrnr_glmnet$new(alpha = 0.50, nfolds = 3)
enet_lrnr_reg75 <- Lrnr_glmnet$new(alpha = 0.75, nfolds = 3)
ranger_lrnr_base <- Lrnr_ranger$new()
ranger_lrnr_ntrees50 <- Lrnr_ranger$new(num.trees = 50)
ranger_lrnr_ntrees100 <- Lrnr_ranger$new(num.trees = 100)
ranger_lrnr_ntrees500 <- Lrnr_ranger$new(num.trees = 500)
xgb_lrnr_base <- Lrnr_xgboost$new(family = "binomial")
xgb50_lrnr <- Lrnr_xgboost$new(nrounds = 50, family = "binomial")
xgb100_lrnr <- Lrnr_xgboost$new(nrounds = 100, family = "binomial")
xgb300_lrnr <- Lrnr_xgboost$new(nrounds = 300, family = "binomial")
dbarts_lrnr_ntree50 <- Lrnr_dbarts$new(family = "binomial", ntree = 50)
dbarts_lrnr_ntree200 <- Lrnr_dbarts$new(family = "binomial", ntree = 200)
dbarts_lrnr_ntree500 <- Lrnr_dbarts$new(family = "binomial", ntree = 500)
hal_lrnr_base <- Lrnr_hal9001$new(fit_type = "glmnet",
                                  family = "binomial",
                                  n_folds = 3)
hal_lrnr_custom <- Lrnr_hal9001$new(fit_type = "glmnet",
                                    family = "binomial",
                                    n_folds = 3,
                                    standardize = FALSE,
                                    lambda.min.ratio = 0.0001,
                                    type.measure = "deviance")
earth_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.earth")
polymars_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.polymars")
nnet_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.nnet")
bayesglm_sl_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.bayesglm")

# meta-learner to ensure predicted probabilities do not go outside [0,1]
logistic_metalearner <- make_learner(Lrnr_solnp, metalearner_logistic_binomial,
                                     loss_loglik_binomial)

# create Super Learner
sl_lrn_reg <- Lrnr_sl$new(
  learners = list(
                  #mean_lrnr,
                  #glm_fast_lrnr,
                  ridge_lrnr,
                  lasso_lrnr,
                  enet_lrnr_reg25,
                  enet_lrnr_reg50,
                  enet_lrnr_reg75,
                  ranger_lrnr_base,
                  ranger_lrnr_ntrees50,
                  ranger_lrnr_ntrees100,
                  ranger_lrnr_ntrees500,
                  xgb_lrnr_base,
                  xgb50_lrnr,
                  xgb100_lrnr,
                  xgb300_lrnr,
                  dbarts_lrnr_ntree50,
                  dbarts_lrnr_ntree200,
                  dbarts_lrnr_ntree500,
                  hal_lrnr_base,
                  hal_lrnr_custom,
                  earth_sl_lrnr,
                  polymars_sl_lrnr,
                  nnet_sl_lrnr,
                  bayesglm_sl_lrnr),
  metalearner = logistic_metalearner
)


# build task for CD4+ outcome regression
cd4_natural_data <- w_baseline %>%
  as_tibble() %>%
  mutate(
    cd4_scaled = janes_data$cd4_score_scaled,
    hiv_status = janes_data$HIVwk28preunbl,
    wts = janes_data$wt_hal
  ) %>%
  as.data.table()

Q_task_cd4_natural <- sl3_Task$new(data = cd4_natural_data,
                                   outcome = "hiv_status",
                                   covariates =
                                     names(cd4_natural_data)[-c(6, 7)],
                                   weights = "wts",
                                   outcome_type = "binomial"
                                  )

# fit SL model for outcome regression in CD4+
sl_train_time_start <- proc.time()
fit_cd4_sl <- sl_lrn_reg$train(Q_task_cd4_natural)
sl_train_time_end <- proc.time()
sl_train_time_cd4 <- sl_train_time_end - sl_train_time_start


# estimate Qn for various shifts
Qn_pred_sl_cd4 <- lapply(delta_shift_grid, function(delta_shift) {
  # create shifted values of intervention
  cd4_shifted_data <- copy(cd4_natural_data) %>%
    as_tibble() %>%
    mutate(
      cd4_scaled = janes_data$cd4_score_scaled + delta_shift
    ) %>%
    as.data.table()
  Q_task_cd4_shifted <- sl3_Task$new(data = cd4_shifted_data,
                                     outcome = "hiv_status",
                                     covariates =
                                       names(cd4_shifted_data)[-c(6, 7)],
                                     weights = "wts",
                                     outcome_type = "binomial"
                                    )
  Qn_natural  <- fit_cd4_sl$predict()
  Qn_shifted  <- fit_cd4_sl$predict(Q_task_cd4_shifted)
  Qn_cd4_sl <- list(
      Qn_natural = Qn_natural,
      Qn_shifted = Qn_shifted
    ) %>%
    as_tibble()

  # output
  return(Qn_cd4_sl)
})
names(Qn_pred_sl_cd4) <- as.character(delta_shift_grid)



# build task for CD8+ outcome regression
cd8_natural_data <- w_baseline %>%
  as_tibble() %>%
  mutate(
    cd8_scaled = janes_data$cd8_score_scaled,
    hiv_status = janes_data$HIVwk28preunbl,
    wts = janes_data$wt_hal
  ) %>%
  as.data.table()

Q_task_cd8_natural <- sl3_Task$new(data = cd8_natural_data,
                                   outcome = "hiv_status",
                                   covariates =
                                     names(cd8_natural_data)[-c(6, 7)],
                                   weights = "wts",
                                   outcome_type = "binomial"
                                  )

# fit SL model for outcome regression in CD8+
sl_train_time_start <- proc.time()
fit_cd8_sl <- sl_lrn_reg$train(Q_task_cd8_natural)
sl_train_time_end <- proc.time()
sl_train_time_cd8 <- sl_train_time_end - sl_train_time_start


# estimate Qn for various shifts
Qn_pred_sl_cd8 <- lapply(delta_shift_grid, function(delta_shift) {
  # create shifted values of intervention
  cd8_shifted_data <- copy(cd8_natural_data) %>%
    as_tibble() %>%
    mutate(
      cd8_scaled = janes_data$cd8_score_scaled + delta_shift
    ) %>%
    as.data.table()
  Q_task_cd8_shifted <- sl3_Task$new(data = cd8_shifted_data,
                                     outcome = "hiv_status",
                                     covariates =
                                       names(cd8_shifted_data)[-c(6, 7)],
                                     weights = "wts",
                                     outcome_type = "binomial"
                                    )
  Qn_natural  <- fit_cd8_sl$predict()
  Qn_shifted  <- fit_cd8_sl$predict(Q_task_cd8_shifted)
  Qn_cd8_sl <- list(
      Qn_natural = Qn_natural,
      Qn_shifted = Qn_shifted
    ) %>%
    as_tibble()

  # output
  return(Qn_cd8_sl)
})
names(Qn_pred_sl_cd8) <- as.character(delta_shift_grid)


# write out nuisance parameter results per intervention node of interest
cd4_nuisance_output <- list(gn = gn_pred_hal_cd4, Qn = Qn_pred_sl_cd4,
                            sl = fit_cd4_sl)
saveRDS(cd4_nuisance_output,
        file = here("data", "cd4_nuisance_delta_grid.rds"))
cd8_nuisance_output <- list(gn = gn_pred_hal_cd8, Qn = Qn_pred_sl_cd8,
                            sl = fit_cd8_sl)
saveRDS(cd8_nuisance_output,
        file = here("data", "cd8_nuisance_delta_grid.rds"))
