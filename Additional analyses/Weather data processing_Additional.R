## Weather data processing (Visual Crossing)
##
## IMPORTANT:
## - The raw daily weather CSVs used by this script were downloaded from the Visual Crossing platform
##   under a time-limited paid subscription.
## - Visual Crossing’s licensing terms indicate that raw data should not be shared/distributed publicly
##   for download; therefore, the raw CSV downloads (Dataframes/weathercsvs/) are not included in this
##   repository. 
## - This script is provided for transparency (i.e., to show exactly how derived, per-site seasonal
##   summaries were computed). The analyses in this repository use those *derived summaries* (and are
##   integrated into the main project dataset), so the main workflow can run without the raw downloads.
##   The below script is provided as an illustrative workflow to show how metrics were computed and derived.

## 0. Setup -------------------------------------------------------------------

# install.packages("here")
library(here)

library(tidyverse)
library(lubridate)
library(zoo)

# List all CSV files
csv_folder <- here::here("Dataframes", "weathercsvs") # Raw data not provided; for illustrative workflow purposes only.
csv_files  <- list.files(csv_folder, pattern = "\\.csv$", full.names = TRUE)

# Initialise an empty list to hold each site’s summary row
summaries <- list()

## 1. Helpers -------------------------------------------------------------------

# Define table to look up each site's respective summer date range ###
site_date_ranges <- tribble(
  ~region,     ~start_date,      ~end_date,
  # Southern hemisphere: 2023‐12 to 2024‐03
  "Argentina", "2023-12-01",     "2024-03-31",
  "Chile",     "2023-12-01",     "2024-03-31",
  "Ecuador",   "2023-12-01",     "2024-03-31",
  "NZ","2023-12-01",    "2024-03-31",
  
  # Northern hemisphere: 2024‐06 to 2024‐09
  "Portugal1", "2024-06-01",     "2024-09-30",
  "Portugal2", "2024-06-01",     "2024-09-30",
  "Vietnam1",  "2024-06-01",     "2024-09-30",
  "Norway",    "2024-06-01",     "2024-09-30",
  "Greenland", "2024-06-01",     "2024-09-30",
  "Denmark",   "2024-06-01",     "2024-09-30",
  
  # UK (summer 2023)
  "UK",        "2023-06-01",     "2023-09-30"
) %>%
  mutate(
    start_date = ymd(start_date),
    end_date   = ymd(end_date)
  )

#  Helper to find runs of ≥3 consecutive days for a given vector
find_runs <- function(dates, logical_vec, min_run = 3) {
  r <- rle(logical_vec)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  tibble(
    start     = dates[starts],
    end       = dates[ends],
    length    = r$lengths,
    condition = r$values
  ) %>%
    filter(condition, length >= min_run) %>%
    mutate(dates = map2(start, end, ~ seq(.x, .y, by = "day")))
}

## 2. Data processing -------------------------------------------------------------------

# Loop over each CSV file
for (file_path in csv_files) {
  # Extract the location name from the filename
  fname   <- basename(file_path)
  region  <- str_split(fname, " - ", simplify = TRUE)[1] %>% str_trim()
  
  # Look up that location’s date range
  date_info <- site_date_ranges %>% filter(region == !!region)
  if (nrow(date_info) == 0) {
    warning(glue::glue("No date‐range found for region '{region}'. Skipping file:\n  {fname}"))
    next
  }
  start_d <- date_info$start_date
  end_d   <- date_info$end_date
  
  # Read the CSV and filter to the correct date range
  df <- read_csv(file_path, show_col_types = FALSE) %>%
    mutate(
      # CSV has column "datetime" that can be parsed as POSIX or character
      date = as_date(datetime)
    ) %>%
    filter(date >= start_d & date <= end_d)
  
  # If after filtering there are zero rows, skip with a warning
  if (nrow(df) == 0) {
    warning(glue::glue("No rows in {region} after date‐filter. Skipping summarization."))
    next
  }
  
  # Compute derived columns
  df <- df %>%
    mutate(
      temp_range    = tempmax - tempmin,
      hw_thresh     = mean(tempmax, na.rm = TRUE) + sd(tempmax, na.rm = TRUE),
      cs_thresh     = mean(tempmin, na.rm = TRUE) - sd(tempmin, na.rm = TRUE),
      heavy_rain_th = quantile(precip, 0.90, na.rm = TRUE),
      gust_th       = quantile(windgust, 0.90, na.rm = TRUE)
    )
  
  # Identify run‐periods for heatwaves, cold snaps, and high UV
  hw_runs <- find_runs(df$date, df$tempmax > df$hw_thresh)
  cs_runs <- find_runs(df$date, df$tempmin < df$cs_thresh)
  uv_runs <- find_runs(df$date, df$uvindex >= 8)
  
  # Summarize all metrics into a tibble
  summary_tbl <- df %>%
    summarize(
      region                 = region,    
      
      # 1. Seasonal means of raw daily variables
      "1_mean_daily_max_temp"    = mean(tempmax,           na.rm = TRUE),
      "1_mean_daily_min_temp"    = mean(tempmin,           na.rm = TRUE),
      "1_mean_daily_temp"        = mean(temp,              na.rm = TRUE),
      "1_mean_humidity"          = mean(humidity,          na.rm = TRUE),
      "1_mean_daily_precip"      = mean(precip,            na.rm = TRUE),
      "1_mean_windgust"          = mean(windgust,          na.rm = TRUE),
      "1_mean_daily_windspeed"   = mean(windspeed,         na.rm = TRUE),
      "1_mean_sealevelpressure"  = mean(sealevelpressure,  na.rm = TRUE),
      "1_mean_cloudcover"        = mean(cloudcover,        na.rm = TRUE),
      "1_mean_visibility"        = mean(visibility,        na.rm = TRUE),
      "1_mean_solarradiation"    = mean(solarradiation,    na.rm = TRUE),
      "1_mean_solarenergy"       = mean(solarenergy,       na.rm = TRUE),
      "1_mean_uvindex"           = mean(uvindex,           na.rm = TRUE),
      
      # 2. Derived variability metrics
      "2_mean_temp_range"        = mean(temp_range,        na.rm = TRUE),
      "2_sd_temp_range"          = sd(temp_range,          na.rm = TRUE),
      "2_sd_precip"              = sd(precip,              na.rm = TRUE),
      "2_sd_windspeed"           = sd(windspeed,           na.rm = TRUE),
      "2_sd_sealevelpressure"    = sd(sealevelpressure,    na.rm = TRUE),
      
      # 3. Derived cumulative metrics
      "3_total_precip"           = sum(precip,             na.rm = TRUE),
      "3_total_solar_energy"     = sum(solarenergy,        na.rm = TRUE),
      
      # 4. Extreme‐event counts & runs
      "4_heatwave_runs"          = nrow(hw_runs),
      "4_heatwave_days"          = sum(hw_runs$length),
      "4_cold_snap_runs"         = nrow(cs_runs),
      "4_cold_snap_days"         = sum(cs_runs$length),
      "4_heavy_rain_days"        = sum(precip > heavy_rain_th,  na.rm = TRUE),
      "4_dry_spell_max"          = max(rle(precip == 0)$lengths),
      "4_gust_extreme_days"      = sum(windgust > gust_th,      na.rm = TRUE),
      "4_sunny_days"             = sum(cloudcover < 20,         na.rm = TRUE),
      "4_high_uv_runs"           = nrow(uv_runs),
      "4_high_uv_days"           = sum(uvindex >= 8,            na.rm = TRUE),
      "4_fog_days"              = sum(visibility < 1,          na.rm = TRUE)
    )
  
  # Append to summaries list
  summaries[[region]] <- summary_tbl
}

# Bind all site summaries into a single tibble (one row per site)
all_sites_summary <- bind_rows(summaries)

view(all_sites_summary)