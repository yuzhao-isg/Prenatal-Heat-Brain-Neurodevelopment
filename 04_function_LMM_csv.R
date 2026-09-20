#-----------------------------------------------------------------------------#
#                         Mixed Model Function -csv                           #
#-----------------------------------------------------------------------------#


###***************************************************************************
###*
###*
###* R code using the Generalized Linear Mixed model to analyze prenatal exposures and outcomes
###* This script also for log transformed outcomes
###* 
###* Created by Yu Zhao, 5/2/2025
###* 
###* 
###***************************************************************************
#' @param data: your dataset
#' @param exposures: exposure variables or list
#' @param outcomes: outcome variables or list
#' @param covariates： covariates
#' @param random_effects: random effect varible/s
#' @param family: using which family, based on the type of your outcomes (e.g., "gaussian", "binomial", "poisson")
#' @param log_transformed_outcome: if the outcomes has been log transformed
#' @param output_dir: save path

run_mixed_models <- function(data, 
                             exposures, 
                             outcomes, 
                             covariates = NULL,
                             random_effects = "(1|id)",
                             family = "gaussian",
                             log_transformed_outcome = FALSE,
                             output_dir = "results/") {
  
  # Load packages
  require(lme4)
  require(dplyr)
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  all_results <- list()
  
  # Single model function
  run_single_model <- function(outcome, exposure, data) {
    
    # Build formula
    formula_str <- paste0(outcome, " ~ ", exposure)
    if (length(covariates) > 0) {
      formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
    }
    formula_str <- paste0(formula_str, " + ", random_effects)
    formula_obj <- as.formula(formula_str)
    
    # Fit model
    model <- tryCatch({
      if (family == "gaussian") {
        lme4::lmer(formula_obj, data = data)
      } else {
        lme4::glmer(formula_obj, data = data, family = family)
      }
    }, error = function(e) NULL)
    
    if (is.null(model)) return(NULL)
    
    # Extract results
    coef_summary <- summary(model)$coefficients
    if (!exposure %in% rownames(coef_summary)) return(NULL)
    
    coef <- coef_summary[exposure, "Estimate"]
    se <- coef_summary[exposure, "Std. Error"]
    
    # Calculate p-value
    if (family == "gaussian") {
      t_val <- coef_summary[exposure, "t value"]
      p_value <- 2 * pt(abs(t_val), df = df.residual(model), lower.tail = FALSE)
    } else {
      p_value <- coef_summary[exposure, "Pr(>|z|)"]
    }
    
    # Calculate CI on log scale
    ci_lower_log <- coef - 1.96 * se
    ci_upper_log <- coef + 1.96 * se
    
    # Keep estimates as coefficients (no transformation)
    if (log_transformed_outcome) {
      # Convert to percent change: 100 × [exp(β) - 1]
      estimate <- 100 * (exp(coef) - 1)
      ci_lower <- 100 * (exp(ci_lower_log) - 1)
      ci_upper <- 100 * (exp(ci_upper_log) - 1)
    } else {
      # Keep raw coefficients
      estimate <- coef
      ci_lower <- ci_lower_log
      ci_upper <- ci_upper_log
    }
    
    # Get sample size
    n <- nobs(model)
    
    return(list(
      outcome = outcome,
      exposure = exposure,
      n = n,
      estimate = estimate,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ))
  }
  
  # Loop through exposures
  for (exposure in exposures) {
    
    exposure_results <- list()
    
    # Run models for each outcome
    for (outcome in outcomes) {
      result <- run_single_model(outcome, exposure, data)
      if (!is.null(result)) {
        exposure_results <- c(exposure_results, list(result))
      }
    }
    
    if (length(exposure_results) == 0) next
    
    # Convert to data frame
    results_df <- do.call(rbind, lapply(exposure_results, function(x) {
      data.frame(
        Exposure = x$exposure,
        Outcome = x$outcome,
        N = x$n,
        Estimate = x$estimate,
        Lower_CI = x$ci_lower,
        Upper_CI = x$ci_upper,
        P_value = x$p_value,
        stringsAsFactors = FALSE
      )
    }))
    
    all_results[[exposure]] <- results_df
  }
  
  # Combine and save all results
  combined_results <- do.call(rbind, all_results)
  write.csv(combined_results, paste0(output_dir, "/all_results.csv"), row.names = FALSE)
  
  return(all_results)
}