
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on October 2025
##' @update Last update July 2026

##' @description Script for main analysis


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


#-----------------------------------------------------------------------------#
#                            CSV file Main model                              #
#-----------------------------------------------------------------------------#

#######################   Neurosonography
load("heat_neurosonography_mixed.RData")

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
                 "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "~/results/Main_results/Neurosonography/",
  log_transformed_outcome = TRUE)

################################################################################

#######################   DP3 28m
load("Heat_DP3_28m_mixed.RData")

results <- run_mixed_models(
  data = Heat_DP3_28m_mixed,
  exposures = c("count_HI_exposure_95_home","log2_count_HI_exposure_95_indoor_day",
                "count_heat_exposure_95_home","log2_count_heat_exposure_95_indoor_day"),
  outcomes = c("ptc_global_development_score","ptc_motricity_score","ptc_adaptive_behaviour_score",
               "ptc_socioemotional_score","ptc_congnition_score","ptc_communication_score"),
  covariates = c("sex_0y_c", "Active_smok_any", "ethnicity_c_2cat", "gestage_0y_c_weeks", 
                 "age_12w_m", "Alcohol_any", "educ_level_m_2cat", "ses_income_acu","dp3_age_28m",
                 "endbf_18m", "ns(as.numeric(conception_date), df = 12)"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "~/results/Main_results/DP3_28m/",
  log_transformed_outcome = FALSE)
###############################################################################



#-----------------------------------------------------------------------------#
#                          Forest Plot Main model                             #
#-----------------------------------------------------------------------------#
rm(list = ls())

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "ggplot2",      
                      "forestplot")

# Set the directory
setwd("~/results/Main_results/")

#######################   Neurosonography

results_neurosonography <- read.csv("Neurosonography/all_results.csv")

# Rename exposure and outcome variables
# Modify this mapping according to your variable names
exposure_names <- c(
  "count_HI_exposure_95_home" = "Outdoor HI 95th",
  "log2_count_HI_exposure_95_indoor_day" = "Indoor HI 95th",
  "count_heat_exposure_95_home" = "Outdoor T 95th",
  "log2_heat_HI_exposure_95_indoor_day" = "Indoor T 95th"
)

neurosonography_names <- c(
  "log_insula" = "Insula depth",
  "log_sylvian_fissure" = "Sylvian fissure depth",
  "log_parieto_occipital_sulcus" = "Parieto-occipital sulcus depth",
  "log_cingulate_sulcus" = "Cingulate sulcus depth",
  "log_calcarine_sulcus" = "Calcarine sulcus depth",
  "log_anterior_ventricle" = "Anterior lateral ventricles width",
  "log_posterior_ventricle" = "Posterior lateral ventricles width",
  "log_third_ventricle" = "Third ventricle width",
  "log_cisterna_magna" = "Cisterna magna width",
  "log_corpus_callosum" = "Corpus callosum length",
  "log_vermis" = "Cerebellar vermis height",
  "log_trans_cerebellar_d" = "Transverse cerebellar diameter"
)

# Add category for outcomes
# Modify this mapping according to your outcomes
neurosonography_categories <- data.frame(
  Outcome_label = c("Insula depth", "Sylvian fissure depth", "Parieto-occipital sulcus depth",
                    "Cingulate sulcus depth", "Calcarine sulcus depth", "Anterior lateral ventricles width",
                    "Posterior lateral ventricles width", "Third ventricle width", "Cisterna magna width",
                    "Corpus callosum length", "Cerebellar vermis height", "Transverse cerebellar diameter"),
  Category = c("Cortical folding", "Cortical folding", "Cortical folding", "Cortical folding","Cortical folding",
               "CSF spaces", "CSF spaces","CSF spaces","CSF spaces",
               "Others","Others","Others"),
  stringsAsFactors = FALSE
)

# Process data
neurosonography_plot <- results_neurosonography %>%
  # Rename exposures
  mutate(Exposure_label = recode(Exposure, !!!exposure_names)) %>%
  # Rename outcomes
  mutate(Outcome_label = recode(Outcome, !!!neurosonography_names)) %>%
  # Add categories
  left_join(neurosonography_categories, by = "Outcome_label") 

category_order <- c("Cortical folding", "CSF spaces", "Others")
outcome_order <- c("Insula depth", "Sylvian fissure depth", "Parieto-occipital sulcus depth",
                   "Cingulate sulcus depth", "Calcarine sulcus depth", "Anterior lateral ventricles width",
                   "Posterior lateral ventricles width", "Third ventricle width", "Cisterna magna width",
                   "Corpus callosum length", "Cerebellar vermis height", "Transverse cerebellar diameter")
neurosonography_plot <- neurosonography_plot %>%
  mutate(Outcome_label = factor(Outcome_label, levels = outcome_order),
         Category = factor(Category, levels = category_order))

# Loop through each exposure to create forest plots
for (exp in unique(neurosonography_plot$Exposure)) {
  
  plot_data <- neurosonography_plot %>%
    filter(Exposure == exp) %>%
    arrange(Category, Outcome_label)
  plot_data$Category <- as.character(plot_data$Category)
  
  exp_label <- unique(plot_data$Exposure_label)
  
  # Create forest dataframe with category grouping
  forest_df <- plot_data %>%
    group_by(Category) %>%
    mutate(first_in_category = row_number() == 1) %>%
    ungroup() %>%
    mutate(
      Category = ifelse(first_in_category, Category, ""),
      Outcome = Outcome_label,
      N = N,
      `  ` = paste(rep(" ", 20), collapse = ""),
      `%Change (95% CI)` = sprintf("%.1f (%.1f, %.1f)", Estimate, Lower_CI, Upper_CI), # keep one decimal
      P = sprintf("%.3f", P_value)
    ) %>%
    select(Category, Outcome, N, `  `, `%Change (95% CI)`, P,
           mean = Estimate, lower = Lower_CI, upper = Upper_CI, first_in_category)
  
  # Create theme with better colors
  tm <- forest_theme(
    base_size = 10,
    core = list(
      bg_params = list(fill = c("#FFFFFF", "#F8F9FA"))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", col = "#2C3E50"),
      bg_params = list(fill = "#D4E6F1")
    )
  )
  
  # Set axis
  xlim <- c(floor(min(forest_df$lower)) - 1, ceiling(max(forest_df$upper)) + 1)
  
  # Create forest plot
  p <- forest(
    forest_df[, 1:6],
    est = forest_df$mean,
    lower = forest_df$lower,
    upper = forest_df$upper,
    sizes = 0.4,
    ci_column = 4,
    xlim = xlim,
    ref_line = 0,
    theme = tm
  )
  
  # Highlight category rows
  for (i in which(forest_df$first_in_category)) {
    p <- edit_plot(p, 
                   row = i + 1,  # +1 for header
                   col = 1,
                   gp = gpar(fontface = "bold", col = "#3498DB", cex = 1.1))
  }
  
  # Save
  filename <- paste0("Neurosonography/forest_P_", gsub("[^A-Za-z0-9]", "_", exp), "1decimal.png")
  png(filename, width = 2800, height = max(1600, nrow(forest_df) * 100), res = 400)
  print(p)
  dev.off()
  
  print(paste("Saved:", filename))
}
################################################################################


#######################    DP3

# Load results
DP3_results <- read.csv("DP3_28m/all_results.csv")

# Rename exposure and outcome variables
# Modify this mapping according to your variable names
exposure_names <- c(
  "count_HI_exposure_95_home" = "Outdoor HI 95th",
  "log2_count_HI_exposure_95_indoor_day" = "Indoor HI 95th",
  "count_heat_exposure_95_home" = "Outdoor T 95th",
  "log2_heat_HI_exposure_95_indoor_day" = "Indoor T 95th"
)

dp3_names <- c(
  "ptc_global_development_score" = "Global development score",
  "ptc_adaptive_behaviour_score" = "Adaptive behaviour score",
  "ptc_motricity_score" = "Motricity score",
  "ptc_socioemotional_score" = "Socioemotional score",
  "ptc_communication_score" = "Communication score",
  "ptc_congnition_score" = "Cognition score"
)

# Process data
dp3_plot <- DP3_results %>%
  # Rename exposures
  mutate(Exposure_label = recode(Exposure, !!!exposure_names)) %>%
  # Rename outcomes
  mutate(Outcome_label = recode(Outcome, !!!dp3_names))

# Loop through each exposure to create forest plots
for (exp in unique(dp3_plot$Exposure)) {
  
  plot_data <- dp3_plot %>%
    filter(Exposure == exp)
  
  exp_label <- unique(plot_data$Exposure_label)
  
  forest_df <- plot_data %>%
    mutate(
      Outcome = Outcome_label,
      N = N,
      `  ` = paste(rep(" ", 20), collapse = ""),
      `Estimate (95% CI)` = sprintf("%.2f (%.2f, %.2f)", Estimate, Lower_CI, Upper_CI), # keep 2 decimal
      P = sprintf("%.3f", P_value)
    ) %>%
    select(Outcome, N, `  `, `Estimate (95% CI)`, P,
           mean = Estimate, lower = Lower_CI, upper = Upper_CI)
  
  # Create theme with better colors
  tm <- forest_theme(
    base_size = 10,
    core = list(
      bg_params = list(fill = c("#FFFFFF", "#F8F9FA"))
    ),
    colhead = list(
      fg_params = list(fontface = "bold", col = "#2C3E50"),
      bg_params = list(fill = "#D4E6F1")
    )
  )
  
  # Set axis
  xlim <- c(floor(min(forest_df$lower)) - 1, ceiling(max(forest_df$upper)) + 1)
  
  # Create forest plot
  p <- forest(
    forest_df[, 1:5],
    est = forest_df$mean,
    lower = forest_df$lower,
    upper = forest_df$upper,
    sizes = 0.4,
    ci_column = 3,
    xlim = xlim,
    ref_line = 0,
    theme = tm
  )
  
  # Save
  filename <- paste0("DP3/forest_P_", gsub("[^A-Za-z0-9]", "_", exp), ".png")
  png(filename, width = 2800, height = max(1600, nrow(forest_df) * 100), res = 500)
  print(p)
  dev.off()
  
  print(paste("Saved:", filename))
}
################################################################################



rm(list = ls())
#-----------------------------------------------------------------------------#
#                        Extended Data :  Neurosonography - DP3               #
#-----------------------------------------------------------------------------#
source("~/script/Function/04_function_LMM_csv.R")
setwd("~/db/output/")
# Load data
load("neurosonography_transvaginal.RData")
load("DP3_28m_covariates_imputed.RData")

setwd("~/results/Main_results/")
# prepare the data
DP3_Neurosonography <- DP3_28m_covariates_imputed %>%
  left_join(neurosonography_transvaginal, by = "id_mother")

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
                 "endbf_18m","parity_m_2cat"),
  random_effects = "(1|hosp_recruit_m_12w)",
  family = "gaussian",
  output_dir = "DP3_Neurosonography/",
  log_transformed_outcome = FALSE)
################################################################################
# END
