library(here)
library(tidyverse)
library(data.table)

# load data sets
fong_data <- read_csv(here("data", "bama.m.for_sharing.csv"))
janes_data <- read_csv(here("data",
                            "v505_tcell_correlates_data_for_sharing.csv"))

# data mangling
## Janes et al. data is in long format, with records indicating different immune
## responses for each observational unit...but has same number of unique IDs as
## Fong et al.
length(fong_data$pub_id) == length(unique(fong_data$pub_id))     ## TRUE
length(fong_data$pub_id) == length(unique(janes_data$pub_id))    ## TRUE
length(fong_data$pub_id) == length(janes_data$pub_id)            ## FALSE


## cellular immune response variables from Janes et al.
## NOTE: should filter Janes et al. data, then put into wide format so as to
##       merge with data from Fong et al.

### first, clean the CD4+ data
janes_data_cd4 <- janes_data %>%
  dplyr::filter(antigen == "ANY VRC ENV" & tcellsub == "CD4+" &
                scoretype == "PolyfunctionalityScore") %>%
  select(-antigen, -scoretype) %>%
  mutate(
    cd4_score = score,
    cd4_score_scaled = score_scaled
  ) %>%
  select(-score, -score_scaled)

### next, the CD8+ data
janes_data_cd8 <- janes_data %>%
  dplyr::filter(antigen == "ANY VRC ENV" & tcellsub == "CD8+" &
                scoretype == "PolyfunctionalityScore") %>%
  select(-antigen, -scoretype) %>%
  transmute(
    pub_id = pub_id,
    cd8_score = score,
    cd8_score_scaled = score_scaled
  )

### now, merge both by subject ID
janes_data_reduced <- left_join(janes_data_cd4, janes_data_cd8, by = "pub_id")


## 6 primary humoral immune responses from Fong et al.
## 1) IgA_env - Summary score for gp120/140 binding
## 2) IgG_env - Summary score for gp120/140 binding
## 3) IgG_V2 - Summary score for V2 binding
## 4) IgG_V3 - Summary score for V3 binding
## 5) IgG_gp41 - gp41 (Subtype B)
## 6) IgG_BioC4_427B - C4 gp120 Peptide (Subtype B)
fong_data_reduced <- fong_data %>%
  select(pub_id, IgA_env, IgG_env, IgG_V2, IgG_V3, IgG_gp41, IgG_BioC4_427B)


## finally, combine Janes et al. and Fong et al. data for our analysis
hvtn505_data <- left_join(janes_data_reduced, fong_data_reduced, by = "pub_id")
saveRDS(hvtn505_data, file = here("data", "hvtn505_for_analysis.rds"))
write_csv(hvtn505_data, path = here("data", "hvtn505_for_analysis.csv"))

