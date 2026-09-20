
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on December 2025
##' @update Last update July 2026

##' @description Script for other extended tables and figures

# Load R libraries ----
list.of.packages <- c("dplyr",
                      "tidyverse",      
                      "tableone",
                      "readr",
                      "corrplot",
                      "patchwork",
                      "ggplot2",
                      "ggmap",
                      "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require,character.only = T) 

# Set the directory
setwd("~/db/")

#---- INDEX ----------------------------------------------------------------                      
# 1)  Supplementary Table 1: Description of socioeconomic, demographic, and lifestyle characteristics of the total BiSC cohort (N = 1080) and the study sample with neurosonography (N = 758) and 28-month follow-up (N = 440).
# 2)  Supplementary Fig. 1: Spearman’s correlation coefficients of the number of heat days across different heat metrics. a. Until neurosonography; b. Until delivery 
# 3)  Extended Data Fig. 2. Histogram of maternal heat exposure days during pregnancy. a. Until neurosonography; b. Until delivery. 
# 4)  Extended Data Fig. 4. Spatial distribution of prenatal outdoor and indoor heat exposure among study participants. a. Until neurosonography; b. Until delivery. 
#---------------------------------------------------------------------------



#-----------------------------------------------------------------------------#
#                         1) Supplementary Table 1                            #
#-----------------------------------------------------------------------------#

# load data
load("output/covariates_neuro_complete.RData")
load("output/neurosonography_covariates_transvaginal_complete.RData")
load("output/DP3_28m_covariates_complete.RData")


bisc_covariates <- covariates_neuro_complete %>%
  select(id_mother, age_12w_m, sex_0y_c, Active_smok_any, gestage_0y_c_weeks, Alcohol_any,
         educ_level_m_2cat, ethnicity_c_2cat) %>%
  mutate(follow = "bisc")

neurosonography_covariates <- neurosonography_covariates_transvaginal_complete %>%
  select(id_mother, age_12w_m, sex_0y_c, Active_smok_any, gestage_0y_c_weeks, Alcohol_any,
         educ_level_m_2cat, ethnicity_c_2cat) %>%
  mutate(follow = "neurosonography") %>%
  filter(!is.na(sex_0y_c))

DP3_28m_covariates <- DP3_28m_covariates_complete  %>%
  select(id_mother, age_12w_m, sex_0y_c, Active_smok_any, gestage_0y_c_weeks, 
         Alcohol_any,educ_level_m_2cat, ethnicity_c_2cat) %>%
  mutate(follow = "DP3_28m")

covariates_neurosonography <- rbind(bisc_covariates, neurosonography_covariates)
covariates_DP3 <- rbind(bisc_covariates, DP3_28m_covariates)

################################################################################

dichotomous_vars <- c("sex_0y_c", "Active_smok_any", "Alcohol_any", "educ_level_m_2cat", "ethnicity_c_2cat")
continuous_vars <- c("age_12w_m", "gestage_0y_c_weeks")
all_vars <- c(dichotomous_vars, continuous_vars)
# Create descriptive table by trimester
table_by_follow_neurosonography <- CreateTableOne(vars = all_vars, strata = "follow", data = covariates_neurosonography, 
                                                  factorVars = dichotomous_vars, includeNA = FALSE)
# Convert to data frame
df_summary_neurosonography <- print(table_by_follow_neurosonography, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")

######
table_by_follow_DP3 <- CreateTableOne(vars = all_vars, strata = "follow", data = covariates_DP3, 
                                      factorVars = dichotomous_vars, includeNA = FALSE)

# Convert to data frame
df_summary_DP3 <- print(table_by_follow_DP3, quote = FALSE, noSpaces = TRUE, printToggle = FALSE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")
###############################################################################


#-----------------------------------------------------------------------------#
#                         2) Supplementary Firgure 1                          #
#-----------------------------------------------------------------------------#
# load data
load("output/exposure/Ultrasound/heat_count_ultrasoundsound.RData")
load("output/exposure/Whole_pregnancy/heat_count_whole_pregnancy.RData")

selected_vars <- c("count_heat_exposure_95_home", "count_HI_exposure_95_home",
                   "count_heat_exposure_95_indoor_day","count_HI_exposure_95_indoor_day")

#######################
### 1) WHOLE PREGNANCY ###
#######################

correlation_selected <- Heat_count_whole[, selected_vars] %>%
  rename(OT95th = count_heat_exposure_95_home,
         OHI95th = count_HI_exposure_95_home,
         IT95th = count_heat_exposure_95_indoor_day,
         IHI95th = count_HI_exposure_95_indoor_day)

correlation_selected <- correlation_selected[, sapply(correlation_selected, is.numeric)]
correlation_selected <- cor(correlation_selected, use = "pairwise.complete.obs")

range(correlation_selected, na.rm = TRUE)
png("~/results/Descriptive/correlation_heat_whole_pregnancy_main.png",
    width = 1200, height = 1000, res = 150)  # Larger size

corrplot(correlation_selected,
         method = "color",
         type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         tl.cex = 1.2,
         number.cex = 1.2)  # Adjust text size if needed

dev.off()

#######################
### 2) ULTRASOUND DATE ###
#######################

correlation_selected <- heat_count_ultrasoundsound[, selected_vars] %>%
  rename(OT95th = count_heat_exposure_95_home,
         OHI95th = count_HI_exposure_95_home,
         IT95th = count_heat_exposure_95_indoor_day,
         IHI95th = count_HI_exposure_95_indoor_day)

correlation_selected <- correlation_selected[, sapply(correlation_selected, is.numeric)]
correlation_selected <- cor(correlation_selected, use = "pairwise.complete.obs")

range(correlation_selected, na.rm = TRUE)
png("~results/Descriptive/correlation_heat_ultrasound_main.png",
    width = 1200, height = 1000, res = 150)  # Larger size

corrplot(correlation_selected,
         method = "color",
         type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         tl.cex = 1.2,
         number.cex = 1.2)  # Adjust text size if needed

dev.off()

#################################################################################


#-----------------------------------------------------------------------------#
#                           3) Extended Data Fig. 2                           #
#-----------------------------------------------------------------------------#

#######################
### 1) WHOLE PREGNANCY ###
#######################

p1_out <- ggplot(Heat_count_whole, aes(x = count_heat_exposure_95_home)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "OT95th", y = "Frequency",
       title = "OT95th") +
  theme_classic()

p2_out <- ggplot(Heat_count_whole, aes(x = count_HI_exposure_95_home)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "OHI95th", y = "Frequency",
       title = "OHI95th") +
  theme_classic()

p1_in <- ggplot(Heat_count_whole, aes(x = count_heat_exposure_95_indoor_day)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "IT95th", y = "Frequency",
       title = "IT95th") +
  theme_classic()

p2_in <- ggplot(Heat_count_whole, aes(x = count_HI_exposure_95_indoor_day)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "IHI95th", y = "Frequency",
       title = "IHI95th") +
  theme_classic()


p_Delivery <- p3_out + p4_out + p3_in + p4_in + plot_layout(ncol = 2, widths = c(1,1)) 

print(p_Delivery)
ggsave("~results/Descriptive/histogram_delivery.png", plot = p_Delivery, width = 10, height = 6, dpi = 600)


#######################
### 2) ULTRASOUND DATE ###
#######################

p3_out <- ggplot(heat_count_ultrasound, aes(x = count_heat_exposure_95_home)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "OT95th", y = "Frequency",
       title = "OT95th") +
  theme_classic()

p4_out <- ggplot(heat_count_ultrasound, aes(x = count_HI_exposure_95_home)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "OHI95th", y = "Frequency",
       title = "OHI95th") +
  theme_classic()

p3_in <- ggplot(heat_count_ultrasound, aes(x = count_heat_exposure_95_indoor_day)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "IT95th", y = "Frequency",
       title = "IT95th") +
  theme_classic()

p4_in <- ggplot(heat_count_ultrasound, aes(x = count_HI_exposure_95_indoor_day)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  labs(x = "IHI95th", y = "Frequency",
       title = "IHI95th") +
  theme_classic()


p_Neurosonography <- p1_out + p2_out + p1_in + p2_in + plot_layout(ncol = 2, widths = c(1,1)) 

print(p_Neurosonography)
ggsave("~results/Descriptive/histogram_neurosonography.png", plot = p_Neurosonography, width = 10, height = 6, dpi = 600)


#-----------------------------------------------------------------------------#
#                           4) Extended Data Fig. 4                           #
#-----------------------------------------------------------------------------#
# load participants coords
load("input/participant_coords.RData")
names(participant)[2] <- "id_mother"

# prepare the map
register_stadiamaps(key = "xxxxxxxxxxxxxxxxxxxxx")
barcelona_bbox <- c(left = 2.05, bottom = 41.30, right = 2.25, top = 41.45)
barcelona_map <- get_stadiamap(bbox = barcelona_bbox, zoom = 12, maptype = "stamen_toner_lite")

#######################
### 1) WHOLE PREGNANCY ###
#######################

load("output/Heat_DP3_28m_mixed.RData")

Heat_DP3_28m_mixed <- Heat_DP3_28m_mixed %>%
  mutate(count_HI_exposure_95_home = count_HI_exposure_95_home * 10,
         count_HI_exposure_95_indoor_day = count_HI_exposure_95_indoor_day * 10,
         count_heat_exposure_95_home = count_heat_exposure_95_home * 10,
         count_heat_exposure_95_indoor_day = count_heat_exposure_95_indoor_day * 10)

participant_heat_day <- Heat_DP3_28m_mixed %>% 
  left_join(participant, by = c("id_mother")) %>%
  select(id_mother,count_HI_exposure_95_home, count_HI_exposure_95_indoor_day, 
         count_heat_exposure_95_home, count_heat_exposure_95_indoor_day,
         longitude, latitude) %>%
  pivot_longer(
    cols = c(count_HI_exposure_95_home, count_HI_exposure_95_indoor_day,
             count_heat_exposure_95_home, count_heat_exposure_95_indoor_day),
    names_to = "variable",
    values_to = "value"
  )

ggmap(barcelona_map) +
  geom_point(
    data = participant_heat_day, 
    aes(x = longitude, y = latitude, color = value), 
    size = 3, alpha = 0.8, shape = 16
  ) +
  scale_color_gradientn(
    colors = c(
      "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"
    ),
    name = "Heat Days Exposure Count",
    trans = "log1p",
    guide = guide_colorbar(barwidth = 15, barheight = 1)
  ) +
  facet_wrap(~ variable, ncol = 2, 
             labeller = as_labeller(c(
               "count_heat_exposure_95_home" = "OT95th",
               "count_heat_exposure_95_indoor_day" = "IT95th",
               "count_HI_exposure_95_home" = "OHI95th",
               "count_HI_exposure_95_indoor_day" = "IHI95th"
             ))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 18)
  )

ggsave("~results/Descriptive/heat_day_map_delivery.png", 
       width = 14, height = 10, dpi = 600) 


#######################
### 2) ULTRASOUND DATE ###
#######################

load("output/heat_neurosonography_mixed.RData")

Heat_Neurosonography_mixed <- Heat_Neurosonography_mixed %>%
  mutate(count_HI_exposure_95_home = count_HI_exposure_95_home * 10,
         count_HI_exposure_95_indoor_day = count_HI_exposure_95_indoor_day * 10,
         count_heat_exposure_95_home = count_heat_exposure_95_home * 10,
         count_heat_exposure_95_indoor_day = count_heat_exposure_95_indoor_day * 10)

participant_heat_day <- Heat_Neurosonography_mixed %>% 
  left_join(participant, by = c("id_mother")) %>%
  select(id_mother,count_HI_exposure_95_home, count_HI_exposure_95_indoor_day, 
         count_heat_exposure_95_home, count_heat_exposure_95_indoor_day,
         longitude, latitude) %>%
  pivot_longer(
    cols = c(count_HI_exposure_95_home, count_HI_exposure_95_indoor_day,
             count_heat_exposure_95_home, count_heat_exposure_95_indoor_day),
    names_to = "variable",
    values_to = "value"
  )

ggmap(barcelona_map) +
  geom_point(
    data = participant_heat_day, 
    aes(x = longitude, y = latitude, color = value), 
    size = 3, alpha = 0.8, shape = 16
  ) +
  scale_color_gradientn(
    colors = c(
      "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"
    ),
    name = "Heat Days Exposure Count",
    trans = "log1p",
    guide = guide_colorbar(barwidth = 15, barheight = 1)
  ) +
  facet_wrap(~ variable, ncol = 2, 
             labeller = as_labeller(c(
               "count_heat_exposure_95_home" = "OT95th",
               "count_heat_exposure_95_indoor_day" = "IT95th",
               "count_HI_exposure_95_home" = "OHI95th",
               "count_HI_exposure_95_indoor_day" = "IHI95th"
             ))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 18)
  )

ggsave("~/results/Descriptive/heat_day_map_neurosonography.png", 
       width = 14, height = 10, dpi = 600) 

################################################################################
# END
