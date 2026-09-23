
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on November 2025
##' @update Last update July 2026

##' @description Script for sensitivity analysis

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "lme4",      
                      "splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

# Set the directory
setwd("~/db/output/")
source("~/script/Function/04_function_LMM_csv.R")

# load data
load("heat_neurosonography_mixed.RData") # Neurosonography
load("Heat_DP3_28m_mixed.RData") # DP3

# load data for complete samples
setwd("~/db/")
load("output/DP3_28m_covariates_complete.RData")
load("output/neurosonography_covariates_transvaginal_complete.RData")
load("exposure/output/Whole_Pregnancy/heat_count_whole_pregnancy.RData")
load("exposure/output/Ultrasound/heat_count_ultrasound.RData")

# load postnatal heat exposure data
load("output/Postnatal/DP3_heat_day_count.RData") 

# load data for IPW analysis
load("output/covariates_neuro_complete.RData")
load("output/neurosonography_transvaginal")
load("output/DP3_28m.RData")

# load data for substudy 3
load("neurosonography_transvaginal.RData")
load("DP3_28m_covariates_imputed.RData")

#---- INDEX ----------------------------------------------------------------
# 01) Complete Samples
#     Supplementary Table 2. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy using all complete cases.
#     Supplementary Table 17. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy using all complete cases.

# 02) Remove Participants with Prenatal Complication
#     Supplementary Table 3. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy using the samples without prenatal complications and preterm birth. 
#     Supplementary Table 18. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy using the samples without prenatal complications and preterm birth.

# 03) IPW Analysis
#      Supplementary Table 4. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy using inverse probability of censoring weights to account for loss to follow-up. 
#      Supplementary Table 19. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy using inverse probability of censoring weights to account for loss to follow-up.

# 04) Alternative Heat Definition: 
#     Supplementary Table 5. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy using alternative definitions. a. Heat index–based exposures (90th and 99th percentiles); b. Temperature-based exposures (90th and 99th percentiles)
#     Supplementary Table 20. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy using alternative definitions. a. Heat index–based exposures (90th and 99th percentiles); b. Temperature-based exposures (90th and 99th percentiles)

# 05) Further Adjust Passive Smoking
#     Supplementary Table 11. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for passive smoking during pregnancy. 
#     Supplementary Table 26. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for passive smoking during pregnancy. 

# 06) Further Adjust Parity
#     Supplementary Table 12. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for parity. 
#     Supplementary Table 27. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for parity. 

# 07) Further Adjust Air Pollution
#     Supplementary Table 14. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for prenatal air pollution (i.e., fine particulate matter and nitrogen dioxide).
#     Supplementary Table 29. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for prenatal air pollution (i.e., fine particulate matter and nitrogen dioxide).

# 08) Further Adjust BDP - ONLY FOR Neurosonography  
#     Supplementary Table 15. Adjusted percent differences in brain morphological structures associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for biparietal diameter.

# 09) Further Adjust Maternal Cognition - ONLY FOR DP3  
#     Supplementary Table 30. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting maternal cognitive performance.

# 10) Further Adjust Postnatal Heat  - ONLY FOR DP3    
#     Supplementary Table 32. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy further adjusting for postnatal heat exposure (number of heat days until DP-3 assessment).

# 11) Excluding Gestational Age 
#     Supplementary Table 31. Adjusteda difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy excluding gestational age as a covariate.

# 12) External_Standardize Score - ONLY FOR DP3
#      Supplementary Table 33. Adjusted difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat-exposure days during pregnancy using the external-standardized scores instead of cohort-specific standardized scores.

# 13) Further Adjust Passive Smoking,  Active Smoking,  Alcohol Consumption - ONLY FOR Sub-study 3
#      Supplementary Table 35. The association between fetal brain morphology and postnatal neurodevelopment was further adjusted for passive smoking during pregnancy*. a. Cortical folding depth; b. Cerebrospinal fluid spaces width; c. Other
#      Supplementary Table 36. The association between fetal brain morphology and postnatal neurodevelopment was further adjusted for active smoking during pregnancy*. a. Cortical folding depth; b. Cerebrospinal fluid spaces width; c. Other
#      Supplementary Table 37. The association between fetal brain morphology and postnatal neurodevelopment was further adjusted for alcohol consumption during pregnancy*. a. Cortical folding depth; b. Cerebrospinal fluid spaces width; c. Other
#---------------------------------------------------------------------------

setwd("~/results/Sensitivity_results/")

#-----------------------------------------------------------------------------#
#                                 01_Complete Samples                         #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

# Prepare the datasets
Heat_Neurosonography_mixed_complete <- heat_count_ultrasound %>%
  mutate(across(
    .cols = -c(id_mother, "log2_count_heat_exposure_90_indoor_day", 
               "log2_count_heat_exposure_95_indoor_day", "log2_count_heat_exposure_99_indoor_day",
               "log2_count_HI_exposure_90_indoor_day","log2_count_HI_exposure_95_indoor_day",       
               "log2_count_HI_exposure_99_indoor_day"), 
    .fns = ~ .x / 10
  )) %>%
  right_join(neurosonography_covariates_transvaginal_complete, by = "id_mother") %>%
  mutate(conception_date = as.Date(conception_date))

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed_complete,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/Neurosonography/01_complete_samples/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

# Prepare the datasets
Heat_DP3_mixed_complete <- heat_count_whole %>%
  mutate(across(
    .cols = -c(id_mother, "log2_count_heat_exposure_90_indoor_day", 
               "log2_count_heat_exposure_95_indoor_day", "log2_count_heat_exposure_99_indoor_day",
               "log2_count_HI_exposure_90_indoor_day","log2_count_HI_exposure_95_indoor_day",       
               "log2_count_HI_exposure_99_indoor_day"), 
    .fns = ~ .x / 10
  )) %>%
  right_join(DP3_28m_covariates_complete, by = "id_mother") %>%
  mutate(conception_date = as.Date(conception_date))

results <- run_mixed_models(
  data = Heat_DP3_mixed_complete,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/01_complete_samples/",
  log_transformed_outcome = FALSE)

################################################################################


#-----------------------------------------------------------------------------#
#                02_Remove Participants with Prenatal Complication            #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

# create a new variable prenatal_complication_combined
Heat_Neurosonography_mixed <- Heat_Neurosonography_mixed %>%
  mutate(Preterm_birth = ifelse(gestage_0y_c_weeks <= 36.99, "yes", "no")) %>%
  mutate(prenatal_complication_combined = ifelse(Ges_Diabetes == "yes" |
                                                   IUGR == "Yes" |
                                                   Hypertention == "yes" |
                                                   Preeclampsia == "yes" |
                                                   Preterm_birth == "yes", "yes", "no"))
without_complications_participants_Neurosonography <- Heat_Neurosonography_mixed %>%
  filter(prenatal_complication_combined == "no")

results <- run_mixed_models(
  data = without_complications_participants_Neurosonography,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/02_remove_prenatal_complication/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

Heat_DP3_28m_mixed <- Heat_DP3_28m_mixed %>%
  mutate(Preterm_birth = ifelse(gestage_0y_c_weeks <= 36.99, "yes", "no")) %>%
  mutate(prenatal_complication_combined = ifelse(Ges_Diabetes == "yes" |
                                                   IUGR == "Yes" |
                                                   Hypertention == "yes" |
                                                   Preeclampsia == "yes" |
                                                   Preterm_birth == "yes", "yes", "no"))
without_complications_participants_DP3 <- Heat_DP3_28m_mixed %>%
  filter(prenatal_complication_combined == "no")

results <- run_mixed_models(
  data = without_complications_participants_DP3,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/02_remove_prenatal_complication/",
  log_transformed_outcome = FALSE)

###############################################################################


#-----------------------------------------------------------------------------#
#                                03_IPW Analysis                              #
#-----------------------------------------------------------------------------#

##### FUNCTION 
#' @param data: your dataset
#' @param exposures: exposure variables or list
#' @param outcomes: outcome variables or list
#' @param covariates： covariates
#' @param random_effects: random effect varible/s
#' @param family: using which family, based on the type of your outcomes (e.g., "gaussian", "binomial", "poisson")
#' @param log_transformed_outcome: if the outcomes has been log transformed
#' @param output_dir: save path
#' 
run_mixed_models_csv_IPW <- function(data, 
                                     exposures, 
                                     outcomes, 
                                     covariates = NULL,
                                     random_effects = "(1|id)",
                                     family = "gaussian",
                                     log_transformed_outcome = FALSE,
                                     output_dir = "results/") {
  
  require(lme4)
  require(dplyr)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  all_results <- list()
  
  run_single_model <- function(outcome, exposure, data, log_transformed_outcome) {
    
    formula_str <- paste0(outcome, " ~ ", exposure)
    if (length(covariates) > 0) {
      formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
    }
    formula_str <- paste0(formula_str, " + ", random_effects)
    formula_obj <- as.formula(formula_str)
    
    model <- tryCatch({
      if (family == "gaussian") {
        lme4::lmer(formula_obj, data = data, weights = weight, subset = R_dp3 == 1)
      } else {
        lme4::glmer(formula_obj, data = data, family = family, weights = weight, subset = R_dp3 == 1)
      }
    }, error = function(e) NULL)
    
    if (is.null(model)) return(NULL)
    
    coef_summary <- summary(model)$coefficients
    if (!exposure %in% rownames(coef_summary)) return(NULL)
    
    coef_val <- coef_summary[exposure, "Estimate"]
    se <- coef_summary[exposure, "Std. Error"]
    
    if (family == "gaussian") {
      t_val <- coef_summary[exposure, "t value"]
      p_value <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
    } else {
      p_value <- coef_summary[exposure, "Pr(>|z|)"]
    }
    
    ci_lower_log <- coef_val - 1.96 * se
    ci_upper_log <- coef_val + 1.96 * se
    
    if (isTRUE(log_transformed_outcome)) {
      estimate <- 100 * (exp(coef_val) - 1)
      ci_lower  <- 100 * (exp(ci_lower_log) - 1)
      ci_upper  <- 100 * (exp(ci_upper_log) - 1)
    } else {
      estimate <- coef_val
      ci_lower  <- ci_lower_log
      ci_upper  <- ci_upper_log
    }
    
    n <- nobs(model)
    
    return(list(
      outcome  = outcome,
      exposure = exposure,
      n        = n,
      estimate = estimate,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value  = p_value
    ))
  }
  
  for (exposure in exposures) {
    
    exposure_results <- list()
    
    for (outcome in outcomes) {
      result <- run_single_model(outcome, exposure, data, log_transformed_outcome)
      if (!is.null(result)) {
        exposure_results <- c(exposure_results, list(result))
      }
    }
    
    if (length(exposure_results) == 0) next
    
    results_df <- do.call(rbind, lapply(exposure_results, function(x) {
      data.frame(
        Exposure  = x$exposure,
        Outcome   = x$outcome,
        N         = x$n,
        Estimate  = x$estimate,
        Lower_CI  = x$ci_lower,
        Upper_CI  = x$ci_upper,
        P_value   = x$p_value,
        stringsAsFactors = FALSE
      )
    }))
    
    all_results[[exposure]] <- results_df
  }
  
  combined_results <- do.call(rbind, all_results)
  write.csv(combined_results, paste0(output_dir, "/all_results.csv"), row.names = FALSE)
  
  return(all_results)
}

#######################
### 1) Neurosonography ###
#######################

# PREPARE THE DATASETS 
neurosonography_transvaginal <- neurosonography_transvaginal %>%
  mutate(R_neuro = 1)
covariates_neurosonography_IPW_heat <- covariates_neuro_complete %>%
  left_join(neurosonography_transvaginal, by = "id_mother") %>%
  mutate(R_neuro = ifelse(is.na(R_neuro),0,1)) %>%
  left_join(heat_count, by = "id_mother") %>%
  mutate(count_heat_exposure_95_home = count_heat_exposure_95_home/10,
         count_HI_exposure_95_home = count_HI_exposure_95_home/10)

##########   probability of being followed up
fit_R_neuro <- glm(R_neuro ~ age_12w_m + educ_level_m_2cat + ethnicity_c_2cat +
                     parity_m_2cat + ses_income_acu + count_heat_exposure_95_home,
                   family = binomial(), data = covariates_neurosonography_IPW_heat)

#########    Predicted probability
covariates_neurosonography_IPW_heat$pr_R_neuro <- predict(fit_R_neuro, type = "response")

########    Calculate the overall follow-up probability
p_R_neuro <- mean(covariates_neurosonography_IPW_heat$R_neuro)

########    Stable weights
# Weights are applied only to those with R=1 (who have an outcome); those with R=0 will not be included in the final analysis.
covariates_neurosonography_IPW_heat$weight <- ifelse(covariates_neurosonography_IPW_heat$R_neuro == 1,
                                                     p_R_neuro / covariates_neurosonography_IPW_heat$pr_R_neuro,NA)

#######    Check the weighted balance
summary(covariates_neurosonography_IPW_heat$weight)
hist(covariates_neurosonography_IPW_heat$weight, breaks = 50)
# weight close to 1 -- Great

######   Final Analysis

results <- run_mixed_models_csv_IPW(
  data = covariates_neurosonography_IPW_heat,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/03_IPW/",
  log_transformed_outcome = TRUE)

################################################################################


#######################
### 2)    DP3    ###
#######################

DP3_28m <- DP3_28m %>%
  mutate(R_dp3 = 1)

covariates_dp3_IPW_heat <- covariates_neuro_complete %>%
  left_join(DP3_28m, by = "id_mother") %>%
  mutate(R_dp3 = ifelse(is.na(R_dp3),0,1)) %>%
  left_join(heat_count, by = "id_mother") %>%
  mutate(delivery_date = as.Date(delivery_date),
         dp3_age_28m = as.numeric(dp3_date_28m - delivery_date)) %>%
  mutate(count_heat_exposure_95_home = count_heat_exposure_95_home/10,
         count_HI_exposure_95_home = count_HI_exposure_95_home/10)

##########   probability of being followed up
fit_R_dp3 <- glm(R_dp3 ~ age_12w_m + educ_level_m_2cat + ethnicity_c_2cat +
                   parity_m_2cat + ses_income_acu + count_heat_exposure_95_home,
                 family = binomial(), data = covariates_dp3_IPW_heat)
#########    Predicted probability
covariates_dp3_IPW_heat$pr_R_dp3 <- predict(fit_R_dp3, type = "response")

########    Calculate the overall follow-up probability
p_R_dp3 <- mean(covariates_dp3_IPW_heat$R_dp3)

########    Stable weights
# Weights are applied only to those with R=1 (who have an outcome); those with R=0 will not be included in the final analysis.
covariates_dp3_IPW_heat$weight <- ifelse(covariates_dp3_IPW_heat$R_dp3 == 1,
                                         p_R_dp3 / covariates_dp3_IPW_heat$pr_R_dp3,NA)

#######    Check the weighted balance
summary(covariates_dp3_IPW_heat$weight)
hist(covariates_dp3_IPW_heat$weight, breaks = 50)
# a little bit right skewness 

results <- run_mixed_models_csv_IPW(
  data = covariates_dp3_IPW_heat,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/03_IPW/",
  log_transformed_outcome = FALSE)
################################################################################


#-----------------------------------------------------------------------------#
#                    04 Alternative Heat Definition                           #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_heat_exposure_90_home","count_HI_exposure_90_home",
                "log2_count_heat_exposure_90_indoor_day","log2_count_HI_exposure_90_indoor_day",
                "count_heat_exposure_99_home","count_HI_exposure_99_home",
                "log2_count_heat_exposure_99_indoor_day","log2_count_HI_exposure_99_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/04_heat_definition/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_90_home","count_HI_exposure_90_home",
                "log2_count_heat_exposure_90_indoor_day","log2_count_HI_exposure_90_indoor_day",
                "count_heat_exposure_99_home","count_HI_exposure_99_home",
                "log2_count_heat_exposure_99_indoor_day","log2_count_HI_exposure_99_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/04_heat_definition/",
  log_transformed_outcome = FALSE)

################################################################################


#-----------------------------------------------------------------------------#
#                         05_Further Adjust Passive Smoking                   #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)", "Passive_smok_any"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/05_adjust_passive_smoking/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)","Passive_smok_any"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3_28m/05_adjust_passive_smoking/",
  log_transformed_outcome = FALSE)


################################################################################

#-----------------------------------------------------------------------------#
#                             06_Further Adjust Parity                        #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)", "parity_m_2cat"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/06_adjust_parity/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)","parity_m_2cat"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/06_adjust_parity/",
  log_transformed_outcome = FALSE)

################################################################################

#-----------------------------------------------------------------------------#
#                       07_Further Adjust Air Pollution                       #
#-----------------------------------------------------------------------------#

#######################
### 1) Neurosonography ###
#######################

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)", "pm25_total_hybrid_32w", "no2_total_hybrid_32w"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/07_adjust_ap/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "pm25_total_hybrid_0y","no2_total_hybrid_0y"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/07_adjust_ap/",
  log_transformed_outcome = FALSE)

###############################################################################


#-----------------------------------------------------------------------------#
#               08_Further Adjust BDP - ONLY FOR Neurosonography              #
#-----------------------------------------------------------------------------#

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)","BDP"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/Neurosonography/07_adjust_BDP/",
  log_transformed_outcome = TRUE)

###############################################################################


#-----------------------------------------------------------------------------#
#              09_Further Adjust Maternal Cognition - ONLY FOR DP3            #
#-----------------------------------------------------------------------------#

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)","pma_tcorrect_m_32w"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/08_adjust_maternal_cognition/",
  log_transformed_outcome = FALSE)

###############################################################################


#-----------------------------------------------------------------------------#
#                10_Further Adjust Postnatal Heat  - ONLY FOR DP3             #
#-----------------------------------------------------------------------------#
Heat_DP3_28m_mixed <- Heat_DP3_28m_mixed %>%
  left_join(DP3_heat_day_count, by = "id_mother") 

# outdoor_T
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "count_post_heat_exposure_95_home"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/10_adjust_postnatal_heat/",
  log_transformed_outcome = FALSE)

# outdoor_HI
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_HI_exposure_95_home"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "count_post_HI_exposure_95_home"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/10_adjust_postnatal_heat/",
  log_transformed_outcome = FALSE)

# Indoor_T
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "count_post_heat_exposure_95_indoor_day"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/10_adjust_postnatal_heat/",
  log_transformed_outcome = FALSE)

# Indoor_HI
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("log2_count_HI_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "count_post_HI_exposure_95_indoor_day"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/10_adjust_postnatal_heat/",
  log_transformed_outcome = FALSE)


###############################################################################


#-----------------------------------------------------------------------------#
#                  11_Excluding Gestational Age - ONLY FOR DP3                #
#-----------------------------------------------------------------------------#

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "/DP3/11_exclude_GA/",
  log_transformed_outcome = FALSE)

###############################################################################

#-----------------------------------------------------------------------------#
#               12_External_Standardize Score - ONLY FOR DP3                  #
#-----------------------------------------------------------------------------#

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day"
  ),
  outcomes = c("pt_global_development_score","pt_motricity_score","pt_adaptive_behaviour_score",
               "pt_socioemotional_score","pt_congnition_score","pt_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "gestage_0y_c_weeks", "ethnicity_c_2cat",
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/12_external_standardize/",
  log_transformed_outcome = FALSE)

###############################################################################


#-----------------------------------------------------------------------------#
#                            13_ONLY FOR SUBSTUDY 3                           #
#-----------------------------------------------------------------------------#

# PREPARE THE DATASETS
DP3_Neurosonography <- DP3_28m_covariates_imputed %>%
  left_join(neurosonography_transvaginal, by = "id_mother")

#######################################
### 1)    Further Adjust Passive Smoking    ###
#######################################
results <- run_mixed_models(
  data = DP3_Neurosonography,
  exposures = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
                "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
                "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
                "log_cingulate_sulcus", "log_insula"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "ethnicity_c_2cat", "gestational_weeks", "gestage_0y_c_weeks",
                 "age_12w_m",  "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m","parity_m_2cat", "Passive_smok_any"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3_Neurosonography/01_adjust_passive_smoking/",
  log_transformed_outcome = FALSE)

#######################################
### 2)    Further Adjust Active Smoking    ###
#######################################
results <- run_mixed_models(
  data = DP3_Neurosonography,
  exposures = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
                "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
                "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
                "log_cingulate_sulcus", "log_insula"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "ethnicity_c_2cat", "gestational_weeks", "gestage_0y_c_weeks",
                 "age_12w_m",  "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m","parity_m_2cat", "Active_smok_any"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3_Neurosonography/02_adjust_active_smoking/",
  log_transformed_outcome = FALSE)

#######################################
### 3)    Further Adjust Alcohol Consumption    ###
#######################################
results <- run_mixed_models(
  data = DP3_Neurosonography,
  exposures = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
                "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
                "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
                "log_cingulate_sulcus", "log_insula"
  ),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "ethnicity_c_2cat", "gestational_weeks", "gestage_0y_c_weeks",
                 "age_12w_m",  "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m","parity_m_2cat", "Alcohol_any"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3_Neurosonography/03_adjust_alcohol_consumption/",
  log_transformed_outcome = FALSE)
################################################################################
# END
