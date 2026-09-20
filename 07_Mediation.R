
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @revisionv1
##' @author Yu Zhao
##' @date Created on 14th August 2026
##' @update Last update 8th September 2026

##' @description Script for mediation analysis
##' 

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "lme4",      
                      "splines",
                      "tidyr",
                      "mediation")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 


#-----------------------------------------------------------------------------#
#                            Mediation analysis                               #
#-----------------------------------------------------------------------------#
# the exposure must occur before the mediator
rm(list = ls())
setwd("~/db/")
load("output/DP3_28m_covariates_imputed.RData")
load("output/neurosonography_transvaginal.RData")
load("exposre/output/heat_count_ultrasound.RData") 


#------------------------------------------------------------
# 1. Variables
#------------------------------------------------------------
exposure <- "log2_count_HI_exposure_95_indoor_day"
mediator <- "log_third_ventricle"
outcome  <- "ptc_congnition_score"
hospital <- "hosp_recruit_m_12w"

common_covariates <- c(
  "sex_0y_c",
  "Active_smok_any",
  "ethnicity_c_2cat",
  "age_12w_m",
  "Alcohol_any",
  "educ_level_m_2cat",
  "ses_income_acu"
)

mediator_covariates <- c(
  common_covariates,
  "gestational_weeks"
)

outcome_covariates <- c(
  common_covariates,
  "gestage_0y_c_weeks",
  "dp3_age_28m",
  "endbf_18m"
)

required_vars <- unique(c(
  exposure,
  mediator,
  outcome,
  hospital,
  mediator_covariates,
  outcome_covariates,
  "conception_date"
))

#------------------------------------------------------------
# 2. Prepare analysis data
#------------------------------------------------------------

analysis_data <- heat_count_ultrasound %>%
  filter(trimester == "Whole_pregnancy") %>%
  dplyr::select(
    id_mother,
    all_of(exposure)
  ) %>%
  inner_join(
    DP3_28m_covariates_imputed,
    by = "id_mother"
  ) %>%
  inner_join(
    neurosonography_transvaginal,
    by = "id_mother"
  ) %>%
  mutate(
    conception_date = as.Date(conception_date),
    conception_date_n = as.numeric(conception_date),
    hosp_recruit_m_12w = factor(hosp_recruit_m_12w)
  ) %>%
  dplyr::select(
    all_of(required_vars),
    conception_date_n
  ) %>%
  drop_na()

#------------------------------------------------------------
# 3. Natural spline for conception date
#------------------------------------------------------------

ns_data <- splines::ns(
  analysis_data$conception_date_n,
  df = 12
) %>%
  as.data.frame()

names(ns_data) <- paste0(
  "ns_cd_",
  seq_len(ncol(ns_data))
)

analysis_data <- bind_cols(
  analysis_data,
  ns_data
)

spline_terms <- names(ns_data)

#------------------------------------------------------------
# 4. Build mixed-model formulas
#------------------------------------------------------------

random_intercept <- paste0(
  "(1 | ",
  hospital,
  ")"
)

formula_m <- reformulate(
  termlabels = c(
    exposure,
    mediator_covariates,
    spline_terms,
    random_intercept
  ),
  response = mediator
)

formula_y <- reformulate(
  termlabels = c(
    exposure,
    mediator,
    outcome_covariates,
    spline_terms,
    random_intercept
  ),
  response = outcome
)

# Display formulas for checking
formula_m
formula_y

#------------------------------------------------------------
# 5. Fit mediator and outcome models
#------------------------------------------------------------

model_m <- lme4::lmer(
  formula = formula_m,
  data = analysis_data,
  REML = FALSE,
  na.action = na.fail
)

model_y <- lme4::lmer(
  formula = formula_y,
  data = analysis_data,
  REML = FALSE,
  na.action = na.fail
)

# Confirm that both models use the same sample
stopifnot(nobs(model_m) == nobs(model_y))

# Inspect model results
summary(model_m)
summary(model_y)

# Check whether the random-effects models are singular
lme4::isSingular(model_m)
lme4::isSingular(model_y)

# Check estimated hospital-level variance
lme4::VarCorr(model_m)
lme4::VarCorr(model_y)

#------------------------------------------------------------
# 6. Multilevel mediation analysis
#------------------------------------------------------------

set.seed(2026)

med_out <- mediation::mediate(
  model.m = model_m,
  model.y = model_y,
  treat = exposure,
  mediator = mediator,
  control.value = 0,
  treat.value = 1,
  group.out = hospital,
  sims = 5000,
  boot = FALSE,
  dropobs = FALSE
)

summary(med_out)
###############################################################################
# END