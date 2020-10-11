library(here)
library(tidyverse)
library(knitr)
library(kableExtra)
library(sl3)

# convenience function to make tables
make_sl_summary <- function(sl_lrnr_fit, file_out) {
  sl_lrnr_risks <- sl_lrnr_fit$cv_risk(loss_squared_error)
  sink(here("tables", paste0(file_out, ".tex")))
  sl_lrnr_risks[, c(1, 2, 6, 3, 7)] %>%
    kable(format = "latex",
          #booktabs = TRUE,
          label = "sl_coefs_risks",
          col.names = c("Learner", "Weight", "Min. Fold Risk", "Mean CV-Risk",
                        "Max. Fold Risk"),
          caption = "Weights and risk estimates of learning algorithms") %>%
    kable_styling() %>%
    print()
  sink()
}

# load data
cd4_data  <- readRDS(here("data", "cd4_nuisance_delta_grid.rds"))
cd8_data  <- readRDS(here("data", "cd8_nuisance_delta_grid.rds"))

# extract SL fits and make tables
cd4_data$sl %>% make_sl_summary(file_out = "cd4_sl_summary")
cd8_data$sl %>% make_sl_summary(file_out = "cd8_sl_summary")
