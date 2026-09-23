# Script: 02_exposure_data.R
##' @project Prenatal heat exposure & Neurodevelopment (neurosonography and DP3)
##' @author Yu Zhao
##' @date Created on November 2025
##' @update Last update July 2026

##' @description 
# Purpose: Prepare heat exposure indicators and cumulative heat-day measures, and outdoor heat wave.
# Inputs: Temperature, participant location, weather-station, and threshold data.
# Outputs: Heat-index data, daily heat indicators, cumulative heat-day, and heat wave files.

# Packages ----
required_packages <- c(
  "dplyr", "lubridate", "sf", "tibble", "tidyr", "weathermetrics"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Set the working directory ----
setwd("~/db/exposure/")
  participant_sf <- sf::st_as_sf(participant, coords = c("longitude", "latitude"), crs = 4326)
# Heat-index calculation ----
# Uses the U.S. National Weather Service heat-index algorithm.

# Whole pregnancy ----
load("input/exposure/temp_hr_bisc.rds")

# Add heat index columns
  # Outdoor
temp_hr_bisc$HI_outdoor_home <- heat.index(
  t = temp_hr_bisc$temperature_mean_home,
  rh = temp_hr_bisc$humidity_mean_home,
  temperature.metric = "celsius",
  round = 2
)
  # Indoor
temp_hr_bisc$HI_indoor_home <- heat.index(
  t = temp_hr_bisc$indoor_temp,
  rh = temp_hr_bisc$indoor_hr,
  temperature.metric = "celsius",
  # Cumulative heat days ----
)

temp_hr_HI_bisc <- temp_hr_bisc
  # Log-transform indoor exposure counts.

# Ultrasound date ----
load("~/db/input/0020_yu_third.RData")
ultrasound_date <- df %>%
  dplyr::select(id_bisc, f_eco_31) %>%
  dplyr::rename(id_mother = id_bisc,
                ultrasound_date = f_eco_31)

temp_hr_HI_bisc_ultrasound <- left_join(temp_hr_HI_bisc, ultrasound_date, by = "id_mother")

temp_hr_HI_bisc_ultrasound <- temp_hr_HI_bisc_ultrasound %>% filter(date <= ultrasound_date)

save(temp_hr_HI_bisc_ultrasound, file = "output/Fixed_day/temp_hr_HI_bisc_ultrasound.rds")

# Fixed calendar day ----
# Prepare the thresholds.

temp_hr_HI_bisc <- temp_hr_HI_bisc %>%
  mutate(month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "warm_season",
           TRUE ~ "cool_season"))


# Indoor threshold based on the participant distribution.

threshold_indoor_day <- temp_hr_HI_bisc %>%
  group_by(month_day) %>%
  summarise(
    Temp_in_w_90th_day = quantile(indoor_temp, 0.90, na.rm = TRUE),
    Temp_in_w_95th_day = quantile(indoor_temp, 0.95, na.rm = TRUE),
    Temp_in_w_99th_day = quantile(indoor_temp, 0.99, na.rm = TRUE),
    HI_in_w_90th_day = quantile(HI_indoor_home, 0.90, na.rm = TRUE),
    HI_in_w_95th_day = quantile(HI_indoor_home, 0.95, na.rm = TRUE),
    HI_in_w_99th_day = quantile(HI_indoor_home, 0.99, na.rm = TRUE)
  )

# Outdoor threshold based on local historical wearther stations
load("input/BNC_threshold.RData")
# Source: https://analisi.transparenciacatalunya.cat/Medi-Ambient/Dades-meteorol-giques-de-la-XEMA/nzvn-apee/about_data
# Missing years or values were obtained from the Meteostat framework.

# Assign each participant to the nearest weather station.
station_coordinate <- read.csv("input/xema_stations.csv") %>%
  rename(station = CODI_ESTACIO, 
         latitude = LATITUD,
         longitude = LONGITUD) %>%
  select(station, latitude, longitude) 
station_sf <- sf::st_as_sf(station_coordinate, coords = c("longitude", "latitude"), crs = 4326)

# participants home
load("input/participant_coords.RData")
participant_sf <- sf::st_as_sf(participant, coords = c("longitude", "latitude"), crs = 4326)

participant_nearest_station <- st_join(participant_sf, station_sf, join = st_nearest_feature) %>%
  select(gid, subject_id, station) %>%
  rename(id_mother = subject_id)
participant_nearest_station <- as.data.frame(participant_nearest_station)
# Define heat days ----

add_heat_indicators <- function(data) {
  data %>%
    mutate(
      heat_exposure_90_home = ifelse(temperature_mean_home > Temp_90th, 1, 0),
      heat_exposure_95_home = ifelse(temperature_mean_home > Temp_95th, 1, 0),
      heat_exposure_99_home = ifelse(temperature_mean_home > Temp_99th, 1, 0),
      HI_exposure_90_home = ifelse(HI_outdoor_home > HI_90th, 1, 0),
      HI_exposure_95_home = ifelse(HI_outdoor_home > HI_95th, 1, 0),
      HI_exposure_99_home = ifelse(HI_outdoor_home > HI_99th, 1, 0),
      heat_exposure_90_indoor_day = ifelse(indoor_temp > Temp_in_w_90th_day, 1, 0),
      heat_exposure_95_indoor_day = ifelse(indoor_temp > Temp_in_w_95th_day, 1, 0),
      heat_exposure_99_indoor_day = ifelse(indoor_temp > Temp_in_w_99th_day, 1, 0),
      HI_exposure_90_indoor_day = ifelse(HI_indoor_home > HI_in_w_90th_day, 1, 0),
      HI_exposure_95_indoor_day = ifelse(HI_indoor_home > HI_in_w_95th_day, 1, 0),
      HI_exposure_99_indoor_day = ifelse(HI_indoor_home > HI_in_w_99th_day, 1, 0)
    )
}

# Whole pregnancy.

heat_define <- temp_hr_HI_bisc %>% 
  left_join(participant_nearest_station, by = c("id_mother", "gid")) %>%
  left_join(threshold, by = c("station","month_day")) %>%
  left_join(threshold_indoor_day, by = "month_day") %>%
  select(-geometry.x, -geometry.y)

# Set an indicator to 1 when the daily temperature or heat index exceeds its threshold.
heat_define <- add_heat_indicators(heat_define)

# Heat days occurred only during the warm season
heat_define <- heat_define %>%
  mutate(
    heat_exposure_90_home = ifelse(season_2 == "warm_season", heat_exposure_90_home, 0),
    heat_exposure_95_home = ifelse(season_2 == "warm_season", heat_exposure_95_home, 0),
    heat_exposure_99_home = ifelse(season_2 == "warm_season", heat_exposure_99_home, 0),
    HI_exposure_90_home = ifelse(season_2 == "warm_season", HI_exposure_90_home, 0),
    HI_exposure_95_home = ifelse(season_2 == "warm_season", HI_exposure_95_home, 0),
    HI_exposure_99_home = ifelse(season_2 == "warm_season", HI_exposure_99_home, 0),
    heat_exposure_90_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_90_indoor_day, 0),
    heat_exposure_95_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_95_indoor_day, 0),
    heat_exposure_99_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_99_indoor_day, 0),
    HI_exposure_90_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_90_indoor_day, 0),
    HI_exposure_95_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_95_indoor_day, 0),
    HI_exposure_99_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_99_indoor_day, 0)
  )

save(heat_define, file = "output/Fixed_day/heat_start_whole_pregnancy.RData")


############      Until Ultrasound Date

# create a month_day variable for our dataset
temp_hr_HI_bisc_ultrasound$MonthDay <- format(temp_hr_HI_bisc_ultrasound$date, "%m-%d")

missing_cols <- setdiff(names(heat_define), names(temp_hr_HI_bisc_ultrasound))
heat_define_ultrasound <- temp_hr_HI_bisc_ultrasound %>%
  left_join(
    heat_define %>% select(id_mother, date, all_of(missing_cols)), 
    by = c("id_mother", "date")
  )

save(heat_define_ultrasound, file = "output/Fixed_day/heat_start_ultrasound.RData")

################################################################################


#-----------------------------------------------------------------------------#
#                          Calculate the cumulative heat days                 #
#-----------------------------------------------------------------------------#

#############      Whole Pregnancy

selected_vars <- c("heat_exposure_90_home","heat_exposure_95_home", "heat_exposure_99_home", 
                   "HI_exposure_90_home", "HI_exposure_95_home","HI_exposure_99_home", 
                   "heat_exposure_90_indoor_day","heat_exposure_95_indoor_day","heat_exposure_99_indoor_day",
                   "HI_exposure_90_indoor_day","HI_exposure_95_indoor_day","HI_exposure_99_indoor_day")

heat_count_whole <- heat_define %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop") 


# log transform for indoor exposure
vas <- c("count_heat_exposure_90_indoor_day","count_heat_exposure_95_indoor_day","count_heat_exposure_99_indoor_day",
         "count_HI_exposure_90_indoor_day","count_HI_exposure_95_indoor_day","count_HI_exposure_99_indoor_day")

heat_count_whole <- heat_count_whole %>%
  mutate(across(all_of(vas), ~log2(.x + 1), .names = "log2_{.col}"))

save(heat_count_whole, file = "output/Fixed_day/heat_count_whole_pregnancy.RData")

############      Until ltrasound Date

heat_count_whole_ultrasound <- heat_define_ultrasound %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop")

# log transform for indoor exposure
heat_count_ultrasound <- heat_count_ultrasound %>%
  mutate(across(all_of(vas), ~log2(.x + 1), .names = "log2_{.col}"))

save(heat_count_ultrasound, file = "output/Fixed_day/heat_count_ultrasound.RData")
# Postnatal heat days were assessed using the same methods and 95th-percentile
# threshold as prenatal heat exposure; the final date was the DP3 measurement date.
################################################################################

# Moving 15-day windows ----

# Indoor threshold.

temp_hr_HI_bisc <- temp_hr_HI_bisc %>%
  mutate(date = as.Date(date),
         month_day = format(date, "%m-%d"),
         season_2 = case_when(
           month_day >= "05-15" & month_day <= "10-15" ~ "warm_season",
           TRUE ~ "cool_season")) 
# reference year
calendar <- tibble(
  target_date = seq(
    as.Date("2018-05-15"),
    as.Date("2018-10-15"),
    by = "day"
  )
) %>%
  mutate(
    target_md = format(target_date, "%m-%d")
  )
window_calendar <- calendar %>%
  rowwise() %>%
  mutate(
    window_dates = list(
      seq(
        target_date - days(7),
        target_date + days(7),
        by = "day"
      )
    )
  ) %>%
  ungroup() %>%
  tidyr::unnest(window_dates) %>%
  mutate(
    window_md = format(window_dates, "%m-%d"),
    target_md = format(target_date, "%m-%d")
  )

threshold_data <- window_calendar %>%
  left_join(
    temp_hr_HI_bisc,
    by = c("window_md" = "month_day")
  )

thresholds_indoor <- threshold_data %>%
  group_by(target_md) %>%
  summarise(
    Temp_in_w_90th_day = quantile(indoor_temp, 0.90, na.rm = TRUE),
    Temp_in_w_95th_day = quantile(indoor_temp, 0.95, na.rm = TRUE),
    Temp_in_w_99th_day = quantile(indoor_temp, 0.99, na.rm = TRUE),
    HI_in_w_90th_day = quantile(HI_indoor_home, 0.90, na.rm = TRUE),
    HI_in_w_95th_day = quantile(HI_indoor_home, 0.95, na.rm = TRUE),
    HI_in_w_99th_day = quantile(HI_indoor_home, 0.99, na.rm = TRUE)
  )


load("input/BNC_threshold_15_moving_window.RData")
# Outdoor threshold.
heat_define <- temp_hr_HI_bisc %>%
  left_join(participant_nearest_station, by = c("id_mother", "gid")) %>%
  left_join(threshold_15_moving_window_outdoor, by = c("station","month_day" = "target_md")) %>%
  left_join(thresholds_indoor, by = c("month_day" = "target_md")) %>%
  select(-geometry.x, -geometry.y)

# Define heat days.
heat_define <- add_heat_indicators(heat_define)
# Heat days occurred only during the warm season
heat_define <- heat_define %>%
  mutate(
    heat_exposure_90_home = ifelse(season_2 == "warm_season", heat_exposure_90_home, 0),
    heat_exposure_95_home = ifelse(season_2 == "warm_season", heat_exposure_95_home, 0),
    heat_exposure_99_home = ifelse(season_2 == "warm_season", heat_exposure_99_home, 0),
    HI_exposure_90_home = ifelse(season_2 == "warm_season", HI_exposure_90_home, 0),
    HI_exposure_95_home = ifelse(season_2 == "warm_season", HI_exposure_95_home, 0),
    HI_exposure_99_home = ifelse(season_2 == "warm_season", HI_exposure_99_home, 0),
    heat_exposure_90_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_90_indoor_day, 0),
    heat_exposure_95_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_95_indoor_day, 0),
    heat_exposure_99_indoor_day = ifelse(season_2 == "warm_season", heat_exposure_99_indoor_day, 0),
    HI_exposure_90_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_90_indoor_day, 0),
    HI_exposure_95_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_95_indoor_day, 0),
    HI_exposure_99_indoor_day = ifelse(season_2 == "warm_season", HI_exposure_99_indoor_day, 0)
  )
save(heat_define, file = "output/Moving_Window/heat_start_whole_pregnancy.RData")

# Count heat days.
selected_vars <- c("heat_exposure_90_home","heat_exposure_95_home", "heat_exposure_99_home", 
                   "HI_exposure_90_home", "HI_exposure_95_home","HI_exposure_99_home", 
                   "heat_exposure_90_indoor_day","heat_exposure_95_indoor_day","heat_exposure_99_indoor_day",
                   "HI_exposure_90_indoor_day","HI_exposure_95_indoor_day","HI_exposure_99_indoor_day")

heat_count_whole <- heat_define %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop") 

# Log-transform indoor and outdoor exposure counts.
vas <- c("count_HI_exposure_90_home", "count_HI_exposure_95_home","count_HI_exposure_99_home",
         "count_heat_exposure_90_home","count_heat_exposure_95_home", "count_heat_exposure_99_home", 
         "count_heat_exposure_90_indoor_day","count_heat_exposure_95_indoor_day","count_heat_exposure_99_indoor_day",
         "count_HI_exposure_90_indoor_day","count_HI_exposure_95_indoor_day","count_HI_exposure_99_indoor_day")

heat_count_whole <- heat_count_whole %>%
  mutate(across(all_of(vas), ~log2(.x + 1), .names = "log2_{.col}"))

save(heat_count_whole, file = "output/Moving_Window/heat_count_whole_pregnancy.RData")

# Ultrasound date.
# Create a month-day variable for the ultrasound dataset.
temp_hr_HI_bisc_ultrasound$MonthDay <- format(temp_hr_HI_bisc_ultrasound$date, "%m-%d")

missing_cols <- setdiff(names(heat_define), names(temp_hr_HI_bisc_ultrasound))
heat_define_ultrasound <- temp_hr_HI_bisc_ultrasound %>%
  left_join(
    heat_define %>% select(id_mother, date, all_of(missing_cols)), 
    by = c("id_mother", "date")
  )

save(heat_define_ultrasound, file = "output/Moving_Window/heat_start_ultrasound.RData")

# Count heat days.
selected_vars <- c("heat_exposure_90_home","heat_exposure_95_home", "heat_exposure_99_home", 
                   "HI_exposure_90_home", "HI_exposure_95_home","HI_exposure_99_home", 
                   "heat_exposure_90_indoor_day","heat_exposure_95_indoor_day","heat_exposure_99_indoor_day",
                   "HI_exposure_90_indoor_day","HI_exposure_95_indoor_day","HI_exposure_99_indoor_day")

heat_count_ultrasound <- heat_define_ultrasound %>%
  group_by(id_mother) %>%
  summarise(across(all_of(selected_vars), ~ sum(. == 1, na.rm = TRUE), .names = "count_{.col}"), .groups = "drop") 

# Log-transform indoor and outdoor exposure counts.
vas <- c("count_HI_exposure_90_home", "count_HI_exposure_95_home","count_HI_exposure_99_home",
         "count_heat_exposure_90_home","count_heat_exposure_95_home", "count_heat_exposure_99_home", 
         "count_heat_exposure_90_indoor_day","count_heat_exposure_95_indoor_day","count_heat_exposure_99_indoor_day",
         "count_HI_exposure_90_indoor_day","count_HI_exposure_95_indoor_day","count_HI_exposure_99_indoor_day")

heat_count_ultrasound <- heat_count_ultrasound %>%
  mutate(across(all_of(vas), ~log2(.x + 1), .names = "log2_{.col}"))

save(heat_count_ultrasound, file = "output/Moving_Window/heat_count_ultrasound.RData")
#################################################################################



#-----------------------------------------------------------------------------#
#                               Outdoor Heat wave                             #
#-----------------------------------------------------------------------------#

# Heat-wave exposure: >=2 consecutive heat days
id_var   <- "id_mother"
date_var <- "date"

heat_vars <- c(
  "HI_exposure_95_home",
  "heat_exposure_95_home"
)
# ------------------------------------------------------------
# Function to calculate heat-wave exposure
# ------------------------------------------------------------

calculate_heatwave <- function(data, heat_var) {
  
  all_ids <- data %>%
    distinct(.data[[id_var]])
  
  episodes <- data %>%
    select(all_of(c(id_var, date_var, heat_var))) %>%
    rename(
      ID = all_of(id_var),
      date = all_of(date_var),
      heat = all_of(heat_var)
    ) %>%
    mutate(
      date = as.Date(date),
      heat = as.numeric(heat)
    ) %>%
    arrange(ID, date) %>%
    group_by(ID) %>%
    mutate(
      consecutive = as.numeric(date - lag(date)) == 1,
      new_episode = heat == 1 &
        (lag(heat, default = 0) == 0 | !consecutive),
      episode = cumsum(new_episode)
    ) %>%
    ungroup() %>%
    filter(heat == 1) %>%
    count(ID, episode, name = "episode_days")
  
  result <- episodes %>%
    group_by(ID) %>%
    summarise(
      # 0/1
      heatwave_2 = as.integer(any(heatwave_days >= 2)),
      # number od heat wave
      heatwave_n_2 = sum(heatwave_days >= 2),
      # number of heat wave days
      heatwave_days_2 = sum(heatwave_days[heatwave_days >= 2],na.rm = TRUE),
      .groups = "drop"
    ) %>%
    right_join(all_ids, by = setNames("ID", id_var)) %>%
    mutate(
      HW2 = replace_na(HW2, 0L),
      HW2_days = replace_na(HW2_days, 0L)
    ) %>%
    rename(
      !!paste0(heat_var, "_HW2") := HW2,
      !!paste0(heat_var, "_HW2_n") := HW2_n,
      !!paste0(heat_var, "_HW2_days") := HW2_days
    )
  
  result
}


# ------------------------------------------------------------
# Calculate heat-wave exposure for a dataset
# ------------------------------------------------------------

calculate_all_heatwaves <- function(data) {
  
  Reduce(
    function(x, y) full_join(x, y, by = id_var),
    lapply(heat_vars, function(x) calculate_heatwave(data, x))
  )
}


# ------------------------------------------------------------
# Whole-pregnancy and ultrasound exposure
# ------------------------------------------------------------

heatwave_whole_pregnancy <- calculate_all_heatwaves(
  heat_define
)

heatwave_ultrasound <- calculate_all_heatwaves(
  heat_define_ultrasound
)


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

output_dir <- "output/Heatwave/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save(
  heatwave_whole_pregnancy,
  file = file.path(output_dir, "heatwave_whole_pregnancy.RData")
)

save(
  heatwave_ultrasound,
  file = file.path(output_dir, "heatwave_ultrasound.RData")
)
################################################################################
# End of script.
