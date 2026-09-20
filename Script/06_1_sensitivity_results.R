
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @revisionv1
##' @author Yu Zhao
##' @date Created on 14th August 2026
##' @update Last update 7th September 2026

##' @description Script for sensitivity analysis

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "lme4",      
                      "splines",
                      "tidyr",
                      "purrr",
                      "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

# Set the directory
setwd("~/db/output/")
source("~/script/Function/04_function_LMM_csv.R")
# load outcome data
load("heat_neurosonography_mixed.RData") # Neurosonography
load("Heat_DP3_28m_mixed.RData") # DP3

# load data for moving window analysis
setwd("~/db/exposure")
load("output/Moving_Window/heat_count_whole_pregnancy.RData")
heat_count_whole_mw <- heat_count_whole
load("output/Moving_Window/heat_count_ultrasound.RData")
heat_count_ultrasound_mw <- heat_count_ultrasound

# load data for heat wave
load("output/Heat_Wave/heatwave_whole_pregnancy.RData")
load("output/Heat_Wave/heatwave_ultrasound.RData")

# load main results for multiple test
setwd("~/results/Main_results/")
results_neurosonography <- read.csv("Neurosonography/all_results.csv")
results_dp3 <- read.csv("DP3_28m/all_results.csv")
#---- INDEX ----------------------------------------------------------------      
# 01) Further adjust for cooling system use
#     Supplementary Table 6. Adjusteda percent differences in brain morphological structures associated with the number of indoor and outdoor heat days during pregnancy further adjusting for frequency of cooling system usage during the warm season.  
#     Supplementary Table 22. Adjusteda difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat days during pregnancy further adjusting for frequency of cooling system usage during the warm season.
# 
# 02) Mutual adjust indoor and outdoor heat days: 
#     Supplementary Table 11. Adjusteda percent differences in brain morphological structures associated with the number of indoor and outdoor heat days during pregnancy in mutually adjusted models. 
#     Supplementary Table 30. Adjusteda difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat days during pregnancy in mutually adjusted models.
#
# 03) Indoor heat as a categorical variables - binary 
#     Supplementary Table 12. Adjusteda percent differences in brain morphological structures associated with categorical cumulative indoor heat-day exposure during pregnancy (0 heat days vs. ≥1 heat days).
#     Supplementary Table 31. Adjusteda difference in postnatal neurodevelopment score associated with categorical cumulative indoor heat-day exposure during pregnancy (0 heat days vs. ≥1 heat days).
#
# 04) Indoor temperature threshold 26℃ 
#    Supplementary Table 13. Adjusteda percent differences in brain morphological structures associated with the number of indoor heat days based on temperature during pregnancy defined using a 26°C threshold.
#    Supplementary Table 32. Adjusteda difference in postnatal neurodevelopment score associated with the number of indoor heat days during pregnancy based on temperature defined using a 26°C threshold.
#
# 05) Outdoor Heat Wave
#    Supplementary Table 15. Adjusteda percent differences in brain morphological structures associated with the outdoor heatwave days during pregnancy.
#    Supplementary Table 34. Adjusteda difference in postnatal neurodevelopment score associated with the outdoor heatwave days during pregnancy.
#
# 06) Threshold: Moving window 15day
#    Supplementary Table 16. Adjusteda percent differences in brain morphological structures associated with the number of outdoor and indoor heat days using 15-day moving-window thresholds.
#    Supplementary Table 35. Adjusteda difference in postnatal neurodevelopment score associated with the number of indoor and outdoor heat days during pregnancy using 15-day moving-window thresholds.
#
# 07) Multiple Test
#    Supplementary Table 17. Adjusteda percent differences in brain morphological structures associated with the number of outdoor and indoor heat days, adjusting the p-value for multiple comparisons. (A) Outdoor heat exposure; (B) Indoor heat exposure
#    Supplementary Table 36. Adjusteda difference in postnatal neurodevelopment score associated with the number of outdoor and indoor heat days, adjusting the p-value for multiple comparisons. (A) Outdoor heat exposure; (B) Indoor heat exposure

#---------------------------------------------------------------------------

setwd("~/results/Sensitivity_results/")
#-----------------------------------------------------------------------------#
#                    01 Further adjust for cooling system use                 #
#-----------------------------------------------------------------------------#
load("db/input/questionnaire_meteo_predictors_adjust.RData")
AC_use <- questionnaire_meteo_predictors_adjust %>%
  distinct(subject_id, .keep_all = TRUE) %>%
  dplyr::select(subject_id, cooling_system_day, cooling_system_night)

Heat_Neurosonography_mixed <- Heat_Neurosonography_mixed %>%
  left_join(AC_use, by = c("id_mother" = "subject_id"))
Heat_DP3_28m_mixed <- Heat_DP3_28m_mixed  %>%
  left_join(AC_use, by = c("id_mother" = "subject_id"))

#######################
### 1) Neurosonography ###
######################
results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_HI_exposure_95_home","log2_count_HI_exposure_95_indoor_day",
                "count_heat_exposure_95_home","log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)", "cooling_system_day","cooling_system_night"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Heat_Neurodevelopment/results/NM_Revise/Adjust_AC_use/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_HI_exposure_95_home","log2_count_HI_exposure_95_indoor_day",
                "count_heat_exposure_95_home","log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)", "cooling_system_day","cooling_system_night"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Heat_Neurodevelopment/results/NM_Revise/Adjust_AC_use/",
  log_transformed_outcome = FALSE)
################################################################################


#-----------------------------------------------------------------------------#
#                 02 Mutual adjust indoor and outdoor heat days               #
#-----------------------------------------------------------------------------#
#######################
### 1) Neurosonography ###
#######################
mutual_adjustment <- list(
  count_heat_exposure_95_home =
    "log2_count_heat_exposure_95_indoor_day",
  log2_count_heat_exposure_95_indoor_day =
    "count_heat_exposure_95_home",
  count_HI_exposure_95_home =
    "log2_count_HI_exposure_95_indoor_day",
  log2_count_HI_exposure_95_indoor_day =
    "count_HI_exposure_95_home"
)
exposures <- names(mutual_adjustment)

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = exposures,
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  exposure_specific_covariates = mutual_adjustment,
  random_effects = "(1 | hosp_recruit_m_12w)",
  family = "gaussian",
  log_transformed_outcome = TRUE,
  output_dir = paste0(
    "Neurosonography/RV1_sensitivity/02_Mutual_adjust"
  )
)

#######################
### 2)    DP3    ###
#######################
results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = exposures,
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  exposure_specific_covariates = mutual_adjustment,
  random_effects = "(1 | hosp_recruit_m_12w)",
  family = "gaussian",
  log_transformed_outcome = TRUE,
  output_dir = paste0(
    "DP3/RV1_sensitivity/02_Mutual_adjust"
  )
)
################################################################################


#-----------------------------------------------------------------------------#
#                03 Indoor heat as a categorical variables - binary           #
#-----------------------------------------------------------------------------#
fit_binary_lmm <- function(
    data,
    exposure,
    outcome,
    covariates,
    random_effect
) {
  
  exposure_term <- paste0(exposure, "_cat")
  
  model_formula <- stats::as.formula(
    paste(
      outcome,
      "~",
      paste(
        c(
          exposure_term,
          covariates,
          paste0("(1 | ", random_effect, ")")
        ),
        collapse = " + "
      )
    )
  )
  
  model <- lmerTest::lmer(
    model_formula,
    data = data,
    REML = FALSE,
    na.action = na.omit
  )
  
  broom.mixed::tidy(
    model,
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    filter(
      term == paste0(exposure_term, "1")
    ) %>%
    transmute(
      exposure = .env$exposure,
      outcome = .env$outcome,
      exposure_group = "1",
      estimate,
      std.error,
      conf.low,
      conf.high,
      p.value
    )
}


# ------------------------------------------------------------------------------
# Run all exposure–outcome combinations
# ------------------------------------------------------------------------------

run_binary_analysis <- function(
    data,
    exposures,
    outcomes,
    covariates,
    random_effect,
    output_file
) {
  
  analysis_data <- data %>%
    mutate(
      across(
        all_of(exposures),
        ~ factor(
          case_when(
            .x == 0 ~ "0",
            .x >= 0.1 ~ "1",
            TRUE ~ NA_character_
          ),
          levels = c("0", "1")
        ),
        .names = "{.col}_cat"
      )
    )
  
  results <- expand_grid(
    exposure = exposures,
    outcome = outcomes
  ) %>%
    mutate(
      result = map2(
        exposure,
        outcome,
        ~ fit_binary_lmm(
          data = analysis_data,
          exposure = .x,
          outcome = .y,
          covariates = covariates,
          random_effect = random_effect
        )
      )
    ) %>%
    select(result) %>%
    unnest(result) %>%
    mutate(
      estimate_CI = sprintf(
        "%.3f (%.3f, %.3f)",
        estimate,
        conf.low,
        conf.high
      )
    )
  
  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write_csv(
    results,
    output_file
  )
  
  results
}

# exposure
indoor_exposures <- c(
  "count_heat_exposure_95_indoor_day",
  "count_HI_exposure_95_indoor_day"
)
#random effects
random_effect <- "hosp_recruit_m_12w"

#######################
### 1) Neurosonography ###
#######################
outcomes <- c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
              "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
              "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
              "log_cingulate_sulcus", "log_insula")
covariates <- c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks",
                "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                "ns(as.numeric(conception_date), df = 12)")

neurosonography_results <- run_binary_analysis(
  data = Heat_Neurosonography_mixed,
  exposures = indoor_exposures,
  outcomes = neurosonography_outcomes,
  covariates = neurosonography_covariates,
  random_effect = random_effect,
  output_file = paste0(
    "Neurosonography/RV1_sensitivity/",
    "03_indoor_category/results.csv"
  )
)

#######################
### 2)    DP3    ###
#######################
outcomes <- c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
              "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score")
covariates <- c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                "endbf_18m", "ns(as.numeric(conception_date), df = 12)")

dp3_results <- run_binary_analysis(
  data = Heat_DP3_28m_mixed,
  exposures = indoor_exposures,
  outcomes = dp3_outcomes,
  covariates = dp3_covariates,
  random_effect = random_effect,
  output_file = paste0(
    "DP3/RV1_sensitivity/03_indoor_category/results.csv"
  )
)
################################################################################


#-----------------------------------------------------------------------------#
#                      04 Indoor temperature threshold 26℃                   #
#-----------------------------------------------------------------------------#
# defined the 26℃ threshold
Indoor_threshold_warm <- temp_hr_HI_bisc %>%
  mutate(date = as.Date(date),
         month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "warm_season",
           TRUE ~ "cool_season")) %>%
  mutate(indoor_26_warm = ifelse(season_2 == "warm_season", indoor_26, 0))

Indoor_threshold_warm_ultra <- temp_hr_HI_bisc_ultrasound %>%
  mutate(date = as.Date(date),
         month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "warm_season",
           TRUE ~ "cool_season")) %>%
  mutate(indoor_26_warm = ifelse(season_2 == "warm_season", indoor_26, 0))

selected_vars <- c("indoor_26_warm")
indoor_count_heat_new <- Indoor_threshold_warm %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop") %>%
  mutate(across(
    .cols = -c(id_mother), 
    .fns = ~ .x / 10
  ))
indoor_count_heat_new_ultra <- Indoor_threshold_warm_ultra %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop") %>%
  mutate(across(
    .cols = -c(id_mother), 
    .fns = ~ .x / 10
  ))

#######################
### 1) Neurosonography ###
#######################

Heat_Neurosonography_mixed <- Heat_Neurosonography_mixed %>%
  left_join(indoor_count_heat_new_ultra, by = "id_mother")

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed,
  exposures = c("count_indoor_24","count_indoor_26",
                "count_indoor_24_warm","count_indoor_26_warm",
                "count_indoor_HI_26.7_warm"),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/RV1_sensitivity/04_indoor_26℃",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################
Heat_DP3_28m_mixed <- Heat_DP3_28m_mixed %>%
  left_join(indoor_count_heat_new, by = "id_mother")

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_indoor_24","count_indoor_26",
                "count_indoor_24_warm","count_indoor_26_warm",
                "count_indoor_HI_26.7_warm"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/RV1_sensitivity/04_indoor_26℃",
  log_transformed_outcome = FALSE)
################################################################################



#-----------------------------------------------------------------------------#
#                              05 Outdoor Heat wave                           #
#-----------------------------------------------------------------------------#
#######################
### 1) Neurosonography ###
#######################
Heat_Neurosonography_mixed_heatwave <- Heat_Neurosonography_mixed %>%
  left_join(heatwave_ultrasound, by = "id_mother")

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed_heatwave,
  exposures = c("HI_exposure_95_home_HW2_days", "heat_exposure_95_home_HW2_days"),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/RV1_sensitivity/05_Heatwave/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################
Heat_DP3_28m_mixed_heatwave <- Heat_DP3_28m_mixed %>%
  left_join(heatwave_whole_pregnancy, by = "id_mother")

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed_heatwave,
  exposures = c("HI_exposure_95_home_HW2_days", "heat_exposure_95_home_HW2_days"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/RV1_sensitivity/05_Heatwave/",
  log_transformed_outcome = FALSE)
################################################################################


#-----------------------------------------------------------------------------#
#                      06 Threshold: Moving window 15day                      #
#-----------------------------------------------------------------------------#
#######################
### 1) Neurosonography ###
#######################
Heat_Neurosonography_mixed_mw <- Heat_Neurosonography_mixed %>%
  left_join(heat_count_ultrasound_mw, by = "id_mother")

results <- run_mixed_models(
  data = Heat_Neurosonography_mixed_mw,
  exposures = c("count_HI_exposure_95_home","log2_count_HI_exposure_95_indoor_day",
                "count_heat_exposure_95_home","log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
               "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
               "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
               "log_cingulate_sulcus", "log_insula"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu",
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "Neurosonography/RV1_sensitivity/06_Moving_Window/",
  log_transformed_outcome = TRUE)

#######################
### 2)    DP3    ###
#######################
Heat_DP3_mixed_mw <- Heat_DP3_mixed %>%
  left_join(heat_count_whole_mw, by = "id_mother")

results <- run_mixed_models(
  data = Heat_DP3_mixed_mw,
  exposures = c("count_HI_exposure_95_home","log2_count_HI_exposure_95_indoor_day",
                "count_heat_exposure_95_home","log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3/RV1_sensitivity/06_Moving_Window/",
  log_transformed_outcome = FALSE)
###############################################################################


#-----------------------------------------------------------------------------#
#                               07 Multiple Tests                             #
#-----------------------------------------------------------------------------#
#######################
### 1) Neurosonography ###
#######################
results_neurosonography <- results_neurosonography %>%
  group_by(Exposure) %>%
  mutate(
    p_BH = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  ungroup()
write.csv(results_neurosonography, 
          file = "Neurosonography/RV1_sensitivity/07_BH.csv",
          row.names = FALSE)

#######################
### 2)    DP3    ###
#######################
results_dp3 <- results_dp3 %>%
  group_by(Exposure) %>%
  mutate(
    p_BH = p.adjust(
      P_value,
      method = "BH"
    )
  ) %>%
  ungroup()
write.csv(results_dp3, 
          file = "DP3/RV1_sensitivity/07_BH.csv",
          row.names = FALSE)
################################################################################