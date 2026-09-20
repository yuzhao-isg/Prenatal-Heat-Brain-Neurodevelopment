
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @revisionv1
##' @author Yu Zhao
##' @date Created on 14th August 2026
##' @update Last update 7th September 2026

##' @description Empirical monitoring-period residual-block resampling

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "lme4",      
                      "splines",
                      "tidyr",
                      "purrr",
                      "broom.mixed",
                      "lmerTest",
                      "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

#===============================================================================
# Temperature-based indoor heat exposure
#===============================================================================
setwd("//FS_HPC.isglobal.lan/HPC_BISC_DATA/analyses/BiSC_24/010_Yu_Zhao/")

# load exposure and outcome database
load("temperature_exposure/results/New_Data/Ultrasound/temp_hr_HI_bisc_ultrasound.rds")
load("Heat_Neurodevelopment/db/Mixed_Linear_Model/heat_neurosonography_mixed.RData")
load("temperature_exposure/results/New_Data/Whole_pregnancy/temp_hr_HI_bisc.rds")
load("Heat_Neurodevelopment/db/Mixed_Linear_Model/Heat_DP3_28m_mixed.RData")
###############  Long-term Temperature
perf_DEML_temp_long <- read_csv("Indoor_model/results/Model/Temporal_Validation_weekly_40/Temperature/temp_perf_DEML_weekly.csv") %>%
  mutate(month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "Warm season",
           TRUE ~ "Cool season"
         )) %>%
  filter(season_2 == "Warm season") %>%
  dplyr::select(id_mother,gid,date,temp_deml_test,indoor_temp, month_day)


###############  Long-term humidity
perf_DEML_hr_long <- read_csv("Indoor_model/results/Model/Temporal_Validation_weekly_40/Humidity/hr_perf_DEML_weekly.csv") %>%
  mutate(month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "Warm season",
           TRUE ~ "Cool season"
         )) %>%
  filter(season_2 == "Warm season") %>%
  dplyr::select(id_mother,gid, date,hr_deml_test,indoor_hr, month_day)

validate_data_long <- left_join(perf_DEML_temp_long, perf_DEML_hr_long, by = c("id_mother", "gid", "date","month_day"))
################################################################################

#-------------------------------------------------------------------------------
# 0. Choose analysis
#-------------------------------------------------------------------------------
# Choose for which analysis: postnatal means DP3 analysis; fetal means neurosonography analysis
# analysis_type <- "postnatal"
analysis_type <- "fetal"

M <- 200
MC_SEED <- 2026

#-------------------------------------------------------------------------------
# 1. Warm-season OOS validation residuals
#-------------------------------------------------------------------------------
validation_temp <- validate_data_long %>%
  mutate(
    date = as.Date(date),
    month_day = format(date, "%m-%d"),
    error_temp = indoor_temp - temp_deml_test
  ) %>%
  filter(
    month_day >= "05-15",
    month_day <= "10-15",
    !is.na(temp_deml_test),
    !is.na(indoor_temp),
    !is.na(error_temp)
  ) %>%
  arrange(id_mother, date)

#-------------------------------------------------------------------------------
# 2. Identify actual contiguous monitoring-period blocks
#-------------------------------------------------------------------------------
validation_temp <- validation_temp %>%
  group_by(id_mother) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    new_block = is.na(lag(date)) | as.numeric(date - lag(date)) > 1,
    validation_block = cumsum(new_block)
  ) %>%
  ungroup()

#-------------------------------------------------------------------------------
# 3. Build empirical residual-block library (1-6 days)
#-------------------------------------------------------------------------------
block_library_temp <- validation_temp %>%
  group_by(id_mother, validation_block) %>%
  summarise(
    block_start = min(date),
    block_end = max(date),
    n_days = n(),
    mean_pred_temp = mean(temp_deml_test, na.rm = TRUE),
    temp_residual_block = list(error_temp),
    .groups = "drop"
  ) %>%
  filter(n_days >= 1, n_days <= 6)

print(table(block_library_temp$n_days))

#-------------------------------------------------------------------------------
# 4. Conditional empirical resampling by predicted-temperature tertile
#-------------------------------------------------------------------------------
temp_block_cuts <- quantile(
  block_library_temp$mean_pred_temp,
  probs = c(1/3, 2/3),
  na.rm = TRUE
)

block_library_temp <- block_library_temp %>%
  mutate(
    temp_stratum = case_when(
      mean_pred_temp <= temp_block_cuts[1] ~ "low",
      mean_pred_temp <= temp_block_cuts[2] ~ "middle",
      TRUE ~ "high"
    )
  )

print(table(block_library_temp$temp_stratum))

#-------------------------------------------------------------------------------
# 5. Full-cohort warm-season predicted indoor temperature + health-model setup
#-------------------------------------------------------------------------------
if (analysis_type == "postnatal") {

  daily <- temp_hr_HI_bisc %>%
    mutate(
      date = as.Date(date),
      month_day = format(date, "%m-%d")
    ) %>%
    filter(month_day >= "05-15", month_day <= "10-15") %>%
    rename(temp_pred = indoor_temp)

  outcomes <- c(
    "ptc_global_development_score",
    "ptc_motricity_score",
    "ptc_adaptive_behaviour_score",
    "ptc_socioemotional_score",
    "ptc_congnition_score",
    "ptc_communication_score"
  )

  covariates <- c(
    "sex_0y_c", "Active_smok_any", "ethnicity_c_2cat",
    "gestage_0y_c_weeks", "age_12w_m", "Alcohol_any",
    "educ_level_m_2cat", "ses_income_acu", "dp3_age_28m",
    "endbf_18m", "ns(as.numeric(conception_date), df = 12)"
  )

  health_dat <- Heat_DP3_28m_mixed

} else if (analysis_type == "fetal") {

  daily <- temp_hr_HI_bisc_ultrasound %>%
    mutate(
      date = as.Date(date),
      month_day = format(date, "%m-%d")
    ) %>%
    filter(month_day >= "05-15", month_day <= "10-15") %>%
    rename(temp_pred = indoor_temp)

  outcomes <- c(
    "log_anterior_ventricle", "log_posterior_ventricle",
    "log_trans_cerebellar_d", "log_vermis", "log_cisterna_magna",
    "log_third_ventricle", "log_corpus_callosum",
    "log_parieto_occipital_sulcus", "log_sylvian_fissure",
    "log_calcarine_sulcus", "log_cingulate_sulcus", "log_insula"
  )

  covariates <- c(
    "sex_0y_c", "Active_smok_any", "ethnicity_c_2cat",
    "gestational_weeks", "age_12w_m", "Alcohol_any",
    "educ_level_m_2cat", "ses_income_acu",
    "ns(as.numeric(conception_date), df = 12)"
  )

  health_dat <- Heat_Neurosonography_mixed

} else {
  stop("analysis_type must be either 'postnatal' or 'fetal'.")
}

daily <- daily %>% arrange(id_mother, date)

#-------------------------------------------------------------------------------
# 6. Heat-day exposure reconstruction
#-------------------------------------------------------------------------------
make_heat_exposure <- function(dat, temp_variable) {

  thresholds <- dat %>%
    group_by(month_day) %>%
    summarise(
      threshold95 = quantile(.data[[temp_variable]], 0.95, na.rm = TRUE),
      .groups = "drop"
    )

  dat %>%
    left_join(thresholds, by = "month_day") %>%
    mutate(
      heat_day_mc = as.integer(.data[[temp_variable]] > threshold95)
    ) %>%
    group_by(id_mother) %>%
    summarise(
      heat_days_mc = sum(heat_day_mc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      heat_log2_mc = log2(heat_days_mc + 1)
    )
}

#-------------------------------------------------------------------------------
# 7. Health-model function
#-------------------------------------------------------------------------------
fit_all_outcomes_lmer <- function(
    exposure_dat,
    health_dat,
    outcomes,
    covariates,
    exposure_var = "heat_log2_mc"
) {

  dat_analysis <- health_dat %>%
    left_join(
      exposure_dat %>% select(id_mother, all_of(exposure_var)),
      by = "id_mother"
    )

  results <- lapply(outcomes, function(y) {

    rhs <- paste(
      c(exposure_var, covariates, "(1 | hosp_recruit_m_12w)"),
      collapse = " + "
    )

    fit <- lmerTest::lmer(
      as.formula(paste(y, "~", rhs)),
      data = dat_analysis,
      REML = FALSE
    )

    tt <- broom.mixed::tidy(
      fit,
      effects = "fixed",
      conf.int = TRUE
    ) %>%
      filter(term == exposure_var)

    data.frame(
      outcome = y,
      beta_raw = tt$estimate,
      se_raw = tt$std.error,
      ci_low_raw = tt$conf.low,
      ci_high_raw = tt$conf.high,
      p = tt$p.value,
      n = nobs(fit)
    )
  })

  bind_rows(results)
}

#-------------------------------------------------------------------------------
# 8. Empirical monitoring-period residual-block resampling
#
# For each participant, residual blocks actually observed in OOS validation are
# sampled WITH REPLACEMENT and sequentially propagated across warm-season daily
# predictions. Internal residual order is preserved. No rnorm() is used.
#-------------------------------------------------------------------------------
simulate_empirical_temp_blocks <- function(
    dat,
    block_library,
    temp_block_cuts
) {

  dat_out <- dat %>%
    arrange(id_mother, date) %>%
    mutate(temp_error_mc = NA_real_)

  ids <- unique(dat_out$id_mother)

  for (id in ids) {

    idx_person <- which(dat_out$id_mother == id)
    idx_person <- idx_person[order(dat_out$date[idx_person])]

    n_person <- length(idx_person)
    pos <- 1

    while (pos <= n_person) {

      # Characterize the next up-to-6 target days by mean predicted temperature.
      lookahead_end <- min(pos + 5, n_person)
      lookahead_idx <- idx_person[pos:lookahead_end]

      mean_target_temp <- mean(
        dat_out$temp_pred[lookahead_idx],
        na.rm = TRUE
      )

      if (is.nan(mean_target_temp)) {
        candidates <- block_library
      } else {

        target_stratum <- dplyr::case_when(
          mean_target_temp <= temp_block_cuts[1] ~ "low",
          mean_target_temp <= temp_block_cuts[2] ~ "middle",
          TRUE ~ "high"
        )

        candidates <- block_library %>%
          filter(temp_stratum == target_stratum)

        if (nrow(candidates) == 0) {
          candidates <- block_library
        }
      }

      # Sample one actual validation residual block with replacement.
      donor_row <- sample(seq_len(nrow(candidates)), size = 1)
      error_vec <- candidates$temp_residual_block[[donor_row]]

      n_remaining <- n_person - pos + 1
      n_assign <- min(length(error_vec), n_remaining)

      # If the pregnancy/warm-season sequence ends before the donor block,
      # retain a random contiguous portion of the donor block.
      if (length(error_vec) > n_assign) {
        start_pos <- sample(
          seq_len(length(error_vec) - n_assign + 1),
          size = 1
        )
        error_vec_use <- error_vec[start_pos:(start_pos + n_assign - 1)]
      } else {
        error_vec_use <- error_vec
      }

      assign_idx <- idx_person[pos:(pos + n_assign - 1)]
      dat_out$temp_error_mc[assign_idx] <- error_vec_use

      pos <- pos + n_assign
    }
  }

  dat_out %>%
    mutate(
      temp_sim = temp_pred + temp_error_mc
    )
}

#-------------------------------------------------------------------------------
# 9. One-simulation sanity check
#-------------------------------------------------------------------------------
set.seed(2026)

test_temp <- simulate_empirical_temp_blocks(
  dat = daily,
  block_library = block_library_temp,
  temp_block_cuts = temp_block_cuts
)

cat("\nMissing simulated residuals:\n")
print(sum(is.na(test_temp$temp_error_mc)))

cat("\nValidation residual distribution:\n")
print(summary(validation_temp$error_temp))

cat("\nResampled residual distribution:\n")
print(summary(test_temp$temp_error_mc))

cat("\nResidual SD comparison:\n")
print(c(
  validation_residual_sd = sd(validation_temp$error_temp, na.rm = TRUE),
  simulated_residual_sd = sd(test_temp$temp_error_mc, na.rm = TRUE)
))

#-------------------------------------------------------------------------------
# 10. Reproduce primary temperature-based exposure/model before MC
#-------------------------------------------------------------------------------
primary_exposure_check <- make_heat_exposure(
  daily,
  temp_variable = "temp_pred"
)

primary_all <- fit_all_outcomes_lmer(
  exposure_dat = primary_exposure_check,
  health_dat = health_dat,
  outcomes = outcomes,
  covariates = covariates
)

print(primary_all)

#-------------------------------------------------------------------------------
# 11. Monte Carlo: empirical temperature residual-block resampling
#-------------------------------------------------------------------------------
run_empirical_temp_blocks <- function(
    M = 200,
    seed = 2030
) {

  set.seed(seed)
  results <- vector("list", M)

  for (m in seq_len(M)) {

    sim_daily <- simulate_empirical_temp_blocks(
      dat = daily,
      block_library = block_library_temp,
      temp_block_cuts = temp_block_cuts
    )

    sim_exposure <- make_heat_exposure(
      sim_daily,
      temp_variable = "temp_sim"
    )

    fits <- fit_all_outcomes_lmer(
      exposure_dat = sim_exposure,
      health_dat = health_dat,
      outcomes = outcomes,
      covariates = covariates,
      exposure_var = "heat_log2_mc"
    )

    fits$iteration <- m
    results[[m]] <- fits

    if (m %% 50 == 0) {
      message("Temperature empirical block simulation: ", m, "/", M)
    }
  }

  bind_rows(results)
}

mc_temp <- run_empirical_temp_blocks(
  M = M,
  seed = MC_SEED
)

# Based on the analysis
# save(mc_temp, file = "Heat_Neurodevelopment/results/NM_Revise/Model_Error/residual_block_dp3.RData")
save(mc_temp, file = "Heat_Neurodevelopment/results/NM_Revise/Model_Error/residual_block_neurosonography.RData")
#-------------------------------------------------------------------------------
# 12. Summarize MC results
# Sim. P2.5-P97.5 = percentile range of simulated point estimates, NOT a 95% CI.
#-------------------------------------------------------------------------------
summary_temp <- mc_temp %>%
  left_join(
    primary_all %>% select(outcome, primary_beta = beta_raw),
    by = "outcome"
  ) %>%
  group_by(outcome) %>%
  summarise(
    primary_beta = first(primary_beta),
    median_beta = median(beta_raw, na.rm = TRUE),
    sim_p025 = quantile(beta_raw, 0.025, na.rm = TRUE),
    sim_p975 = quantile(beta_raw, 0.975, na.rm = TRUE),
    same_direction_pct = mean(
      sign(beta_raw) == sign(first(primary_beta)),
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  )

#-------------------------------------------------------------------------------
# 13. Reporting scale
#-------------------------------------------------------------------------------
if (analysis_type == "fetal") {

  empirical_table_temp <- summary_temp %>%
    mutate(
      primary_estimate = 100 * (exp(primary_beta) - 1),
      median_sim_estimate = 100 * (exp(median_beta) - 1),
      sim_p025_estimate = 100 * (exp(sim_p025) - 1),
      sim_p975_estimate = 100 * (exp(sim_p975) - 1)
    ) %>%
    select(
      outcome,
      primary_estimate,
      median_sim_estimate,
      sim_p025_estimate,
      sim_p975_estimate,
      same_direction_pct
    )

} else {

  empirical_table_temp <- summary_temp %>%
    transmute(
      outcome,
      primary_estimate = primary_beta,
      median_sim_estimate = median_beta,
      sim_p025_estimate = sim_p025,
      sim_p975_estimate = sim_p975,
      same_direction_pct
    )
}

print(empirical_table_temp)

#-------------------------------------------------------------------------------
# 14. Optional save
#-------------------------------------------------------------------------------
# write.csv(
#   empirical_table_temp,
#   paste0("temperature_empirical_block_MC_", analysis_type, ".csv"),
#   row.names = FALSE
# )

#===============================================================================
# End
#===============================================================================
