
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on September 2025
##' @update Last update July 2026

##' @description Script for GAM analysis for linearity assessment


# Load R libraries ----
list.of.packages <- c("dplyr",
                      "mgcv",      
                      "splines",
                      "itsadug")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

# Set the directory
setwd("~/db/")

#-----------------------------------------------------------------------------#
#                              Prepare Database                               #
#-----------------------------------------------------------------------------#
load("exposure/output/Fixed_day/heat_count_whole_pregnancy.RData")
load("exposure/output/Fixed_day/heat_count_ultrasound.RData")
load("output/neurosonography_covariates_transvaginal_imputed.RData")
load("output/DP3_28m_covariates_imputed.RData")

###############  Neurosonography
Heat_Neurosonography_mixed <- heat_count_ultrasound %>%
  mutate(across(
    .cols = -c(id_mother, "log2_count_heat_exposure_90_indoor_day", 
               "log2_count_heat_exposure_95_indoor_day", "log2_count_heat_exposure_99_indoor_day",
               "log2_count_HI_exposure_90_indoor_day","log2_count_HI_exposure_95_indoor_day",       
               "log2_count_HI_exposure_99_indoor_day"), 
    .fns = ~ .x / 10 # unit increase ten days
  )) %>%
  right_join(neurosonography_covariates_transvaginal_imputed, by = "id_mother") %>%
  mutate(conception_date = as.Date(conception_date)) 

save(Heat_Neurosonography_mixed, file = "output/heat_neurosonography_mixed.RData")


###############  DP3
Heat_DP3_mixed <- heat_count_whole %>%
  mutate(across(
    .cols = -c(id_mother, "log2_count_heat_exposure_90_indoor_day", 
               "log2_count_heat_exposure_95_indoor_day", "log2_count_heat_exposure_99_indoor_day",
               "log2_count_HI_exposure_90_indoor_day","log2_count_HI_exposure_95_indoor_day",       
               "log2_count_HI_exposure_99_indoor_day"), 
    .fns = ~ .x / 10 # unit increase ten days
  )) %>%
  right_join(DP3_28m_covariates_imputed, by = "id_mother") %>%
  mutate(conception_date = as.Date(conception_date)) 

save(Heat_DP3_mixed, file = "output/heat_dp3_mixed.RData")
###############################################################################


#-----------------------------------------------------------------------------#
#                                 GAM Function                                #
#-----------------------------------------------------------------------------#
##' @param data: dataset
##' @param outcomes: outcome variables or list
##' @param exposures: exposure variables or list
##' @param covariates： covariates list
##' @param output_dir: save path
##' @param k_spline: Basis dimension (k) for the spline smooth term in the GAM.

fit_and_plot_gam <- function(
    data,
    outcomes,
    exposures,
    covariates,                   
    output_dir = "gam_results",
    k_spline = 10                 
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  covariate_string <- paste(covariates, collapse = " + ")
  
  for (outcome in outcomes) {
    for (exposure in exposures) {
      
      message("Fitting model for outcome: ", outcome, ", exposure: ", exposure)
      
      # formula
      formula_str <- paste0(
        outcome, " ~ s(", exposure, ", k=", k_spline, ") + ",
        covariate_string
      )
      
      formula_obj <- as.formula(formula_str)
      
      # gamm
      model <- gamm(
        formula_obj,
        random = list(hosp_recruit_m_12w = ~1),
        data = data
      )
      
      # Save summary
      summary_path <- file.path(output_dir, paste0("summary_", outcome, "_", exposure, ".txt"))
      sink(summary_path)
      print(summary(model$gam))
      cat("\n\n--- GAM Model Check ---\n")
      gam.check(model$gam)
      sink()
      
      # Save plot
      plot_path <- file.path(output_dir, paste0("plot_", outcome, "_", exposure, ".png"))
      png(plot_path, width = 800, height = 600, res = 120)
      plot(
        model$gam, pages = 1, se = TRUE, shade = TRUE,
        xlab = exposure,
        ylab = "Partial effect",
        main = paste(outcome)
      )
      dev.off()

    }
  }
  
  message("All models have been fitted and saved in: ", output_dir)
}

################################################################################


#-----------------------------------------------------------------------------#
#                                    GAM                                      #
#-----------------------------------------------------------------------------#

####################  Neurosonography
covariates_neurosonography <- c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestational_weeks", 
                                "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu", 
                                "ns(as.numeric(conception_date), df = 12)")

# Define exposure variables
exposure_variables <- c("count_heat_exposure_95_home","count_HI_exposure_95_home",
                        "count_heat_exposure_95_indoor_day","count_HI_exposure_95_indoor_day",
                        "log2_count_heat_exposure_95_indoor_day","log2_count_HI_exposure_95_indoor_day")

# Define outcome variables 
outcome_neurosonography <- c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
                       "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
                       "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
                       "log_cingulate_sulcus", "log_insula")

save_dir_neurosonography <- "~/results/GAM/Neurosonography/"

Results <- fit_and_plot_gam(data = Heat_Neurosonography_mixed, outcomes = outcome_neurosonography,
                            exposures = exposure_variables, covariates = covariates_neurosonography,
                            output_dir = save_dir_neurosonography)


####################  DP3 

covariates_dp3 <- c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                        "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                        "endbf_18m", "ns(as.numeric(conception_date), df = 12)")

outcome_dp3 <- c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
                       "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score")

save_dir_dp3 <- "~/results/GAM/DP3/"

Results <- fit_and_plot_gam(data = Heat_DP3_28m_mixed, outcomes = outcome_dp3,
                            exposures = exposure_variables, covariates = covariates_dp3,
                            output_dir = save_dir_dp3)
################################################################################


####################   DP3 and neurosonography
load("output/DP3_28m_covariates_imputed.RData")
load("output/neurosonography_transvaginal.RData")

DP3_28m_Neurosonography <- DP3_28m_covariates_imputed %>%
  left_join(neurosonography_transvaginal, by = "id_mother")

Neurosonography <- c("log_anterior_ventricle", "log_posterior_ventricle", "log_trans_cerebellar_d",
                        "log_vermis", "log_cisterna_magna", "log_third_ventricle", "log_corpus_callosum",
                        "log_parieto_occipital_sulcus", "log_sylvian_fissure", "log_calcarine_sulcus",
                        "log_cingulate_sulcus", "log_insula")

covariates <- c("sex_0y_c", "ethnicity_c_2cat", "gestage_0y_c_weeks", "gestational_weeks",
                        "age_12w_m", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                        "endbf_18m")

dp3 <- c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
                       "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score")

save_dir <- "~/results/GAM/neurosonography_dp3/"
Results <- fit_and_plot_gam(data = DP3_28m_Neurosonography, outcomes = dp3,
                            exposures = Neurosonography, covariates = covariates,
                            output_dir = save_dir)

#################################################################################
# END