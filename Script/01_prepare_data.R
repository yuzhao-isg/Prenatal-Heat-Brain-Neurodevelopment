
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on March 2025
##' @update Last update July 2026

##' @description 
# Purpose: Preparing the outcomes and covariates dataset, including loading data and data imputation
# Inputs: Original variables from BiSC database
# Outputs: Outcomes and covariates in the analysis


# Load R libraries ----
list.of.packages <- c("dplyr",
                      "tidyverse",      
                      "mice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

# Set the directory
setwd("~/db/")

#-----------------------------------------------------------------------------#
#                   Prepare the Neurosonography Datasets                      #
#-----------------------------------------------------------------------------#

# Load neurosonography data 
load("input/2024_10_17_Yu_Zhao_BiSC_24_010.RData")
# load rater information
load("input/2024_11_11_Yu_Examiner.RData")
# load ultrasound date
load("input/0020_yu_third.RData")
ultrasound_date <- df

neurosonography <- yu_bisc_24_010_2024_10_17 %>%
  select(id_mother, id_child, present_cef_32w, VLDa_31, VLIa_31, atrium_c_32w, atrium_laterality_c_32w,
         dtc_31, vermis_31, Cm_31, IIIv_31, dfo_31, cc_31, ic_31, ceph_index_c_32w, cpov_31,
         val_cis_poc, gra_cis_poc, csv_31, val_cis_Sil, gra_cis_sil, ccal_31,
         val_cis_cal, gra_cis_cal, ccin_31, val_cis_cin, gra_cis_cin, ins_31,
         val_ins, cvsa_31, osa_31)

# For the Anterior lateral ventricles width, calculate the average value (left and right)
neurosonography$anterior_ventricle <- rowMeans(neurosonography[,c("VLDa_31","VLIa_31")], na.rm = TRUE)
neurosonography <- neurosonography %>%
  rename(posterior_ventricle = atrium_c_32w,
         trans_cerebellar_d = dtc_31,
         vermis = vermis_31,
         cisterna_magna = Cm_31,
         third_ventricle = IIIv_31,
         corpus_callosum = cc_31,
         parieto_occipital_sulcus = val_cis_poc,
         sylvian_fissure = val_cis_Sil,
         calcarine_sulcus = val_cis_cal,
         cingulate_sulcus = val_cis_cin,
         insula = val_ins,
         right_anterior_v = VLDa_31,
         left_anterior_v = VLIa_31)

log_transform_vars <- c("anterior_ventricle", "posterior_ventricle", "trans_cerebellar_d",
                        "vermis", "cisterna_magna","third_ventricle","corpus_callosum",
                        "parieto_occipital_sulcus", "sylvian_fissure", "calcarine_sulcus",
                        "cingulate_sulcus", "insula")


neurosonography_log <- neurosonography %>%
  mutate(across(
    all_of(log_transform_vars),
    ~ if_else(.x == 0, NA_real_, log(.x)),
    .names = "log_{.col}"
  ))

neurosonography_log_rater <- left_join(neurosonography_log, yu_examiner, by = "id_mother")

# add the variable of ultrasound date
ultrasound_date <- ultrasound_date %>%
  dplyr::rename(id_mother = id_bisc,
                ultrasound_date = f_eco_31,
                BDP = dbp_31,
                gestational_weeks = gestage_ultrasound_weeks_m_32w,
                gestational_days = gestage_ultrasound_m_32w,
                head_cir = PC_31) %>%
  select(id_mother, ultrasound_date)

# Only keep the participants with transvaginal approach
neurosonography_transvaginal <- neurosonography_log_rater %>%
  filter(present_cef_32w == "cefalica") %>%
  left_join(ultrasound_date, by = "id_mother")

save(neurosonography_transvaginal, file = "output/neurosonography_transvaginal.RData")
################################################################################

#-----------------------------------------------------------------------------#
#                          Prepare the DP3 Datasets                           #
#-----------------------------------------------------------------------------#
load("input/DP3_c_28m_v2_AS_20240802.RData")

# clean the datsets, set -999 as NA
dp3_28m[sapply(dp3_28m, is.numeric)] <- lapply(dp3_28m[sapply(dp3_28m, is.numeric)], 
                                               function(x) ifelse(x == -999, NA, x))
DP3_28m <- dp3_28m %>%
  rename(id_mother = id_bisc,
         dp3_date_28m = Fecha) %>%
  mutate(dp3_date_28m = as.Date(dp3_date_28m)) %>%
  rename(pt_inconsistency = pt_Inconsistencia,
         pt_global_development_score =  `pt_Índice general de desarrollo`,
         pt_motricity_score = `pt_Motricidad`,
         pt_adaptive_behaviour_score = `pt_Cond. adaptativa`,
         pt_socioemotional_score = `pt_Socioemocional`,
         pt_congnition_score = `pt_Cognición`,
         pt_communication_score = `pt_Comunicación`,
         pd_inconsistency = pd_Inconsistencia,
         pd_global_development_score =  `pd_Índice general de desarrollo`,
         pd_motricity_score = `pd_Motricidad`,
         pd_adaptive_behaviour_score = `pd_Cond. adaptativa`,
         pd_socioemotional_score = `pd_Socioemocional`,
         pd_congnition_score = `pd_Cognición`,
         pd_communication_score = `pd_Comunicación`) %>%
  # transform to cohort-specific standardized scores
  mutate(across(
    c(pd_motricity_score,pd_adaptive_behaviour_score,
      pd_socioemotional_score,pd_congnition_score,pd_communication_score), 
    ~ ((.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)) * 15 + 100,
    .names = "{str_replace(.col, 'pd_', 'ptc_')}"
  )) %>%
  mutate(ptc_global_development_score = (ptc_motricity_score + ptc_adaptive_behaviour_score +
                                           ptc_socioemotional_score + ptc_congnition_score + ptc_communication_score)/5) 

save(DP3_28m, file = "output/DP3_28m.RData")
################################################################################


#-----------------------------------------------------------------------------#
#                             Prepare Covariates                              #
#-----------------------------------------------------------------------------#
load("input/2024_10_17_Yu_Zhao_BiSC_24_010.RData")
load("input/0014_yuzhao.RData")
covariates_1 <- df
load("input/alcohol_yu.RData") # reported alcohol comsumption during pregnancy
alcohol <- df
load("input/efw_32w.RData") # estimated fetal weight
estimated_fetal_weight <- df %>%
  rename(id_mother = id)
# load other ultrasound information
load("input/0020_yu_third.RData")
ultrasound_infor <- df
# load the breast feeding at 18m
load("input/bayley18m_dataset.RData")
bf <- df
bf$id <- as.numeric(as.character(bf$id))
# load maternal cognitive perform
load("input/pma_tcorrect_m_32w.RData")
Maternal_PMA <- df %>%
  rename(id_mother = id_bisc)

colnames(covariates_1)
colnames(yu_bisc_24_010_2024_10_17)

covariates_1 <- covariates_1 %>%
  select(id_bisc, age_12w_m, educ_level_m_3cat, ethnicity_c_4cat, sex_0y_c,
         parity_m_2cat, BISC_cb_v01_12w_questMareEV12w.smoke_any_m, hosp_recruit_m_12w,
         bmi_m_12w, gestage_0y_c_weeks, weight_0y_c, delivery_type, C_DG, C_hip,
         C_ret, fgr_c, hdp_denovo_m, pe_total_m, covid_frontier_m_3c, hospital,
         BISC_cb_v01_12w_questMareEV12w.smokepassive_m_preg_2c,
         weight_0y_c) %>%
  left_join(alcohol, by = "id_bisc") %>%
  rename(id_mother = id_bisc,
         Ges_Diabetes = C_DG,
         hospital_delivery = hospital,
         Pre_Complic = C_hip,
         IUGR = C_ret,
         FGR = fgr_c,
         Hypertention = hdp_denovo_m,
         Preeclampsia = pe_total_m,
         Active_smok_any = BISC_cb_v01_12w_questMareEV12w.smoke_any_m,
         Passive_smok_any = BISC_cb_v01_12w_questMareEV12w.smokepassive_m_preg_2c,
         Alcohol_any = alcohol_m_preg_2c,
         birth_weight = weight_0y_c) %>%
  left_join(estimated_fetal_weight, by = "id_mother")

covariates_1 <- covariates_1 %>%
  mutate(hospital_delivery_3c = ifelse(hospital_delivery == "HSjD" | hospital_delivery == "Maternitat",
                                       "HSjD_Maternitat", ifelse(hospital_delivery == "Home" |
                                                                   hospital_delivery == "Other", "Others", "HSPau"))) %>%
  mutate(educ_level_m_2cat = ifelse(educ_level_m_3cat == "University", "with_university", "without_university")) %>%
  mutate(ethnicity_c_2cat = ifelse(ethnicity_c_4cat == "eur", "eur", "other")) %>%
  mutate(rater = ifelse(hosp_recruit_m_12w == "MT" | hosp_recruit_m_12w == "SJD", "MT-SJD", "SP"))

# add the variables of air pollution 
covariates_2 <- yu_bisc_24_010_2024_10_17 %>%
  select(id_mother,id_child,ses_income_acu, no2_total_hybrid_32w,
         pm25_total_hybrid_32w, no2_total_hybrid_0y, 
         pm25_total_hybrid_0y,  fur_eco_1, f_parto) %>%
  rename(conception_date = fur_eco_1,
         delivery_date = f_parto) %>%
  mutate(conception_month = month(conception_date))

covariates <- left_join(covariates_1, covariates_2, by = "id_mother")

# Convert some variables to factors
factor_vars <- c("educ_level_m_3cat", "ethnicity_c_4cat", "sex_0y_c", "parity_m_2cat",
                 "Active_smok_any", "hosp_recruit_m_12w", "delivery_type", "Ges_Diabetes",
                 "Pre_Complic", "IUGR", "FGR", "Hypertention", "Preeclampsia", 
                 "Passive_smok_any", "Alcohol_any", "season_conception_m",
                 "hospital_delivery_3c","educ_level_m_2cat", "ethnicity_c_2cat", "rater")
covariates[factor_vars] <- lapply(covariates[factor_vars], as.factor)

# Based on previous email, 10030111 10045911 11001911 12034311 10052711 are Girls 
covariates_complete <- covariates %>%
  mutate(sex_0y_c = as.character(sex_0y_c)) %>%
  mutate(sex_0y_c = ifelse(id_mother %in% c("10030111", "10045911", "11001911", "12034311", "10052711"), "Girl", sex_0y_c))

# breasting feeding
bf <- bf %>%
  select("id", "endbf_18m") %>%
  rename(id_mother = id)

# add other ultrasound information
ultrasound_infor <- ultrasound_infor %>%
  dplyr::rename(id_mother = id_bisc,
                BDP = dbp_31,
                gestational_weeks = gestage_ultrasound_weeks_m_32w,
                gestational_days = gestage_ultrasound_m_32w,
                head_cir = PC_31) %>%
  select(id_mother, gestational_weeks, gestational_days, BDP, head_cir)

covariates_neuro_complete <- covariates_complete %>%
  left_join(bf, by = "id_mother") %>%
  left_join(Maternal_PMA, by = "id_mother") %>%
  left_join(ultrasound_infor, by = "id_mother")

save(covariates_neuro_complete, file = "output/covariates_neuro_complete.RData")
#################################################################################


#-----------------------------------------------------------------------------#
#                             multiple imputation                             #
#-----------------------------------------------------------------------------#

covariates_neuro_complete <- covariates_neuro_complete%>%
  filter(!is.na(sex_0y_c))
# 12037711 with the neurosonography (transabdominal) information but no sex 


##variables used to predict missings (from the other variables)
features <- c("educ_level_m_3cat", "ethnicity_c_4cat", "birth_weight", "hosp_recruit_m_12w",
              "delivery_type",  "season_conception_m","BDP", "head_cir",
              "hospital_delivery_3c","sex_0y_c", "rater", "ethnicity_c_2cat",
              "ses_income_acu","Ges_Diabetes", "Pre_Complic","IUGR", "FGR", 
              "Hypertention", "Preeclampsia","gestational_days",
              "age_12w_m","parity_m_2cat","gestage_0y_c_weeks")

##variables to impute
to_impute <- c("Active_smok_any","bmi_m_12w", "Passive_smok_any", 
               "Alcohol_any", "efw_32w", "endbf_18m", "pma_tcorrect_m_32w")


apply(covariates_neuro_complete[, to_impute], 2, function(x) sum(is.na(x))) ##number of missings
apply(covariates_neuro_complete[, to_impute], 2, function(x) round(100*sum(is.na(x))/nrow(covariates_neuro_complete),2)) 

################################################################################

id <- c("id_mother") # id variable
data <- covariates_neuro_complete[, c(id, features, to_impute)] #subset only to 
data[, id] <- NULL

predictorMatrix <- quickpred(data, mincor = 0.2, minpuc = 0.4)
prova_covariates <- mice(data,  m = 100, print = F, seed = 1234, predictorMatrix = predictorMatrix)
covariates_neuro_imputed <- complete(prova_covariates,include = F)

###evaluate which one database is best
covariates_neuro_imputed <- cbind(covariates_neuro_complete[, id], covariates_neuro_imputed)

covariates_neuro_complete[, to_impute] <- covariates_neuro_imputed[, to_impute] ## replace variables by imputed values
covariates_neuro_complete[, features] <- covariates_neuro_imputed[, features] ## replace variables by features values
covariates_neuro_imputed <- covariates_neuro_complete

save(covariates_neuro_imputed, file = "output/covariates_neuro_imputed.RData")

################################################################################


#-----------------------------------------------------------------------------#
#                  Prepare the Outcome + Covariates Database                  #
#-----------------------------------------------------------------------------#

##############    neurosonography
neurosonography_covariates_transvaginal_imputed <- neurosonography_transvaginal %>%
  left_join(covariates_neuro_imputed, by = "id_mother")

neurosonography_covariates_transvaginal_complete <- neurosonography_transvaginal %>%
  left_join(covariates_neuro_complete, by = "id_mother")

save(neurosonography_covariates_transvaginal_imputed, file = "output/neurosonography_covariates_transvaginal_imputed.RData")
save(neurosonography_covariates_transvaginal_complete, file = "output/neurosonography_covariates_transvaginal_complete.RData")


############# DP3
DP3_28m_covariates_imputed <- DP3_28m %>%
  left_join(covariates_neuro_imputed, by = "id_mother") %>%
  mutate(delivery_date = as.Date(delivery_date),
         dp3_age_28m = as.numeric(dp3_date_28m - delivery_date) / 30.4375)

DP3_28m_covariates_complete <- DP3_28m %>%
  left_join(covariates_neuro_complete, by = "id_mother") %>%
  mutate(delivery_date = as.Date(delivery_date),
         dp3_age_28m = as.numeric(dp3_date_28m - delivery_date) / 30.4375)

save(DP3_28m_covariates_imputed, file = "output/DP3_28m_covariates_imputed.RData")
save(DP3_28m_covariates_complete, file = "output/DP3_28m_covariates_complete.RData")
##################################################################################
# END

