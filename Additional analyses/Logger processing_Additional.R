## 0. Setup -------------------------------------------------------------------

# install.packages("here")   # Install once per machine should code and data be read from GitHub
library(here)

library(tidyverse)
library(lubridate)
library(patchwork)
library(ggtext)
library(grid)
library(zoo)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read and process data function
read_and_process_data <- function(directory, skip_rows, site_label, pattern = "\\.csv$") {
  csv_files <- list.files(
    path = directory,
    recursive = TRUE,
    pattern = pattern,
    full.names = TRUE
  )
  
  data_frames_list <- list()
  for (csv_file in csv_files) {
    df_name <- tools::file_path_sans_ext(basename(csv_file))
    df <- read.csv(csv_file, skip = skip_rows, stringsAsFactors = FALSE)
    if(!"logger_name" %in% names(df)){
      df$logger_name <- df_name
    }
    df <- df %>% mutate(site = site_label)
    data_frames_list[[df_name]] <- df
  }
  
  logger_data_final <- bind_rows(data_frames_list, .id = "source_id")
  
  # Convert "time" to date-time and split into date and time components
  logger_data_final$time <- ymd_hms(logger_data_final$time)
  logger_data_final <- logger_data_final %>%
    separate(time, into = c("date", "clock_time"), sep = " ") %>%
    separate(clock_time, into = c("hour", "min", "sec"), sep = ":")
  
  # Define pad groups based on the logger name
  logger_data_final <- logger_data_final %>%
    mutate(
      log_group = case_when(
        grepl("W", logger_name) ~ "white_logs",
        grepl("B", logger_name) ~ "black_logs",
        TRUE ~ NA_character_
      ),
      Plate_ID = logger_name,
      group = case_when(
        log_group == "white_logs" ~ "Control pads",
        log_group == "black_logs" ~ "Warm pads",
        TRUE ~ NA_character_
      )
    )
  
  return(logger_data_final)
}

## 1. Helpers -------------------------------------------------------------------

# Return the maximum value of x (ignoring NAs) or NA if x is all NA
max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

# For each unique logger (Plate_ID), find overall maximum temperature
# Group by site and pad group to compute the mean of those maximums
summarize_max_metric <- function(data, metric_name) {
  data %>%
    filter(!is.na(Plate_ID)) %>%
    group_by(Plate_ID) %>%
    summarise(
      metric = max_or_na(temp),
      site = first(site),
      group = first(group),
      .groups = "drop"
    ) %>%
    group_by(site, group) %>%
    summarise(!!metric_name := mean(metric, na.rm = TRUE), .groups = "drop")
}

# For each logger, calculate the daily maximum temperature
# Then averages these daily maximums (ADM) and compute a standard error by site and group
summarize_daily_metric_with_se <- function(data, metric_name) {
  df_pad <- data %>%
    filter(!is.na(Plate_ID)) %>%
    group_by(Plate_ID, date) %>%
    summarise(daily_max = max_or_na(temp), .groups = "drop") %>%
    left_join(
      dplyr::distinct(dplyr::select(data, Plate_ID, group, site)),
      by = "Plate_ID"
    ) %>%
    group_by(Plate_ID) %>%
    summarise(
      ADM = mean(daily_max, na.rm = TRUE),
      group = first(group),
      site = first(site),
      .groups = "drop"
    )
  
  df_pad %>%
    group_by(site, group) %>%
    summarise(
      mean_ADM = mean(ADM, na.rm = TRUE),
      se_ADM = sd(ADM, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

## 2. Process each individual location -------------------------------------------------------------------
# Note; some warning messages appear due to trailing rows and summary/footer lines. NAs are dropped later in code.

# Southern hemisphere sites
  data_argentina <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Argentina_loggers"),
    skip_rows = 21,
    site_label = "Argentina"
  )
  
  data_chile <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Chile_loggers"),
    skip_rows = 21,
    site_label = "Chile"
  )
  
  data_ecuador <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Ecuador_loggers"),
    skip_rows = 21,
    site_label = "Ecuador"
  )
  
  data_nz <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "NZ_loggers"),
    skip_rows = 20,
    site_label = "New Zealand"
  )

# Northern hemisphere sites
  data_portugal_1 <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Portugal1_loggers"),
    skip_rows = 21,
    site_label = "Portugal 1",
    pattern = "P.*\\.csv$"  # Only include files with "P"
  )
  
  data_portugal_2 <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Portugal2_loggers"),
    skip_rows = 21,
    site_label = "Portugal 2"
  )
  
  data_denmark <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Denmark_loggers"),
    skip_rows = 21,
    site_label = "Denmark"
  )
  
  data_greenland <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Greenland_loggers"),
    skip_rows = 21,
    site_label = "Greenland"
  )

  data_norway <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Norway_loggers"),
    skip_rows = 21,
    site_label = "Norway"
  )
  
  data_vietnam <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "Vietnam_loggers"),
    skip_rows = 21,
    site_label = "Vietnam"
  )

# UK (summer 2023)
  data_UK <- read_and_process_data(
    directory = here::here("Dataframes", "HotMess loggers", "UK_loggers"),
    skip_rows = 20,
    site_label = "UK"
  )

# Combine sites
southern_data_all <- bind_rows(data_argentina, data_chile, data_ecuador, data_nz)
northern_data_all <- bind_rows(data_portugal_1, data_portugal_2, data_denmark, data_greenland, data_norway, data_vietnam)

# Define summer periods (matching periods for 'Weather data processing' code also):
#   Southern summer: 2023-12-01 to 2024-03-31
#   Northern summer: 2024-06-01 to 2024-09-30
#   UK (northern summer of 2023)
southern_summer <- southern_data_all %>%
  filter(date > as.Date("2023-12-01") & date < as.Date("2024-03-31"))
northern_summer <- northern_data_all %>%
  filter(date > as.Date("2024-06-01") & date < as.Date("2024-09-30"))
UK_summer <- data_UK %>%
  filter(date > as.Date("2023-06-01") & date < as.Date("2023-09-30"))

## 3. Calculate summaries for each metric for each site -------------------------------------------------------------------

# Southern hemisphere sites
summary_summer_max_south <- summarize_max_metric(southern_summer, "Mean_Max_Summer")
summary_summer_adm_south <- summarize_daily_metric_with_se(southern_summer, "ADM_Summer")

# Northern hemisphere sites

summary_summer_max_north <- summarize_max_metric(northern_summer, "Mean_Max_Summer")
summary_summer_adm_north <- summarize_daily_metric_with_se(northern_summer, "ADM_Summer")

# UK (2023)
summary_summer_max_UK <- summarize_max_metric(UK_summer, "Mean_Max_Summer")
summary_summer_adm_UK <- summarize_daily_metric_with_se(UK_summer, "ADM_Summer")

# Recode group names to "Warmed" and "Ambient"
recode_groups <- function(df) {
  df %>%
    mutate(group = case_when(
      group == "Warm pads" ~ "Warmed",
      group == "Control pads" ~ "Ambient",
      TRUE ~ group
    )) %>%
    mutate(group = factor(group, levels = c("Warmed", "Ambient")))
}

summary_summer_max_south <- recode_groups(summary_summer_max_south)
summary_summer_adm_south <- recode_groups(summary_summer_adm_south)
summary_summer_max_north <- recode_groups(summary_summer_max_north)
summary_summer_adm_north <- recode_groups(summary_summer_adm_north)
summary_summer_max_UK <- recode_groups(summary_summer_max_UK)
summary_summer_adm_UK <- recode_groups(summary_summer_adm_UK)

south_summary <- full_join(summary_summer_max_south, summary_summer_adm_south, by = c("site", "group"))
north_summary <- full_join(summary_summer_max_north, summary_summer_adm_north, by = c("site", "group"))
UK_summary <- full_join(summary_summer_max_UK, summary_summer_adm_UK, by = c("site", "group"))
summary_table_all <- bind_rows(south_summary, north_summary, UK_summary) %>%
  rename(
    `Mean max` = Mean_Max_Summer,
    `Mean ADM` = mean_ADM,
    `SE Mean ADM` = se_ADM
  )
print(summary_table_all)

## 4. Plotting for each site -------------------------------------------------------------------

# Create a minimal ggplot object with the site name as a title, used for labeling panels
make_site_title_plot <- function(site_name, font_size = 5) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = site_name, 
             size = font_size, fontface = "bold") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_site_panel <- function(site_name, summary_max, summary_adm) {
  # Filter data for the current site
  summary_max_site <- summary_max %>% filter(site == site_name)
  summary_adm_site <- summary_adm %>% filter(site == site_name)
  
  ## Dynamic y-limits for mean max
  mm_vals <- summary_max_site$Mean_Max_Summer
  if (!all(is.na(mm_vals))) {
    mm_min <- min(mm_vals, na.rm = TRUE)
    mm_max <- max(mm_vals, na.rm = TRUE)
    mm_range <- mm_max - mm_min
    buffer_mm <- 0.05 * mm_range  
    if (mm_range == 0) {
      mm_min <- mm_min - 0.5
      mm_max <- mm_max + 0.5
    } else {
      mm_min <- mm_min - buffer_mm
      mm_max <- mm_max + buffer_mm
    }
    mm_limits <- c(mm_min, mm_max)
  } else {
    mm_limits <- NULL
  }
  
  ## Dynamic y-limits for mean ADM (including error bars)
  adm_vals <- summary_adm_site$mean_ADM
  adm_errs <- summary_adm_site$se_ADM
  if (!all(is.na(adm_vals))) {
    lower_vals <- adm_vals - adm_errs
    upper_vals <- adm_vals + adm_errs
    adm_min <- min(lower_vals, na.rm = TRUE)
    adm_max <- max(upper_vals, na.rm = TRUE)
    adm_range <- adm_max - adm_min
    buffer_adm <- 0.05 * adm_range  
    if (adm_range == 0) {
      adm_min <- adm_min - 0.5
      adm_max <- adm_max + 0.5
    } else {
      adm_min <- adm_min - buffer_adm
      adm_max <- adm_max + buffer_adm
    }
    adm_limits <- c(adm_min, adm_max)
  } else {
    adm_limits <- NULL
  }
  
  # Create mean max bar plot
  p_mean_max <- ggplot(summary_max_site, aes(x = group, y = Mean_Max_Summer, fill = group)) +
    geom_col(width = 0.8) +
    scale_fill_manual(name = "Group", values = c("Warmed" = "orange", "Ambient" = "grey75")) +
    labs(title = "Mean max", x = "Group", y = "Temperature (°C)") +
    (if (!is.null(mm_limits)) coord_cartesian(ylim = mm_limits) else NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "none",
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8))
  
  # Create mean ADM bar plot
  p_adm <- ggplot(summary_adm_site, aes(x = group, y = mean_ADM, fill = group)) +
    geom_col(width = 0.8) +
    geom_errorbar(aes(ymin = mean_ADM - se_ADM, ymax = mean_ADM + se_ADM),
                  width = 0.2, size = 0.5, color = "black") +
    scale_fill_manual(name = "Group", values = c("Warmed" = "orange", "Ambient" = "grey75")) +
    labs(title = "Mean ADM", x = "Group", y = "Temperature (°C)") +
    (if (!is.null(adm_limits)) coord_cartesian(ylim = adm_limits) else NULL) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "none",
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 8))
  
  # Create site title plot
  title_plot <- make_site_title_plot(site_name, font_size = 5)
  
  # Combine the two bar plots
  bars <- p_mean_max + p_adm
  
  # Stack the site title above plots
  site_panel <- wrap_plots(title_plot, bars, ncol = 1, heights = c(0.15, 1))
  
  return(site_panel)
}

# Create individual panels for each site
panel_argentina  <- make_site_panel("Argentina", summary_summer_max_south, summary_summer_adm_south)
panel_chile      <- make_site_panel("Chile",     summary_summer_max_south, summary_summer_adm_south)
panel_ecuador    <- make_site_panel("Ecuador",   summary_summer_max_south, summary_summer_adm_south)
panel_nz         <- make_site_panel("New Zealand", summary_summer_max_south, summary_summer_adm_south)
panel_portugal_1 <- make_site_panel("Portugal 1",  summary_summer_max_north, summary_summer_adm_north)
panel_portugal_2 <- make_site_panel("Portugal 2",  summary_summer_max_north, summary_summer_adm_north)
panel_denmark    <- make_site_panel("Denmark",  summary_summer_max_north, summary_summer_adm_north)
panel_greenland  <- make_site_panel("Greenland",  summary_summer_max_north, summary_summer_adm_north)
panel_norway     <- make_site_panel("Norway",  summary_summer_max_north, summary_summer_adm_north)
panel_UK         <- make_site_panel("UK",  summary_summer_max_UK, summary_summer_adm_UK)
panel_vietnam   <- make_site_panel("Vietnam",  summary_summer_max_north, summary_summer_adm_north)


# Combine the panels into a grid
southern_panel <- wrap_plots(
  list(panel_argentina, panel_ecuador, panel_chile, panel_nz, ggplot() + theme_void(), ggplot() + theme_void()),
  ncol = 2, nrow = 4
)

northern_panel <- wrap_plots(
  list(panel_portugal_1, panel_portugal_2, panel_denmark, panel_greenland, panel_norway, panel_UK, panel_vietnam),
  ncol = 2, nrow = 4
)

final_panel <- northern_panel | southern_panel
final_panel

## 5. Calculate summaries for each metric globally -------------------------------------------------------------------

# Compute global summary across all sites by group
global_summary <- summary_table_all %>%
  group_by(group) %>%
  summarise(
    Global_Mean_max = mean(`Mean max`, na.rm = TRUE),
    Global_Mean_ADM = mean(`Mean ADM`, na.rm = TRUE)
  ) %>%
  ungroup()

# Pivot to wide format so that Warmed and Ambient are separate
global_summary_wide <- global_summary %>%
  pivot_wider(names_from = group, values_from = c(Global_Mean_max, Global_Mean_ADM))

# Calculate overall raw differences
global_raw_overall <- tibble(
  Metric = c("Mean max", "Mean ADM"),
  Raw_Difference = c(
    global_summary_wide$Global_Mean_max_Warmed - global_summary_wide$Global_Mean_max_Ambient,
    global_summary_wide$Global_Mean_ADM_Warmed - global_summary_wide$Global_Mean_ADM_Ambient
  )
)

# Calculate overall proportional differences (exploratory)
global_prop_overall <- tibble(
  Metric = c("Mean max", "Mean ADM"),
  Proportional_Difference = c(
    (global_summary_wide$Global_Mean_max_Warmed - global_summary_wide$Global_Mean_max_Ambient) / global_summary_wide$Global_Mean_max_Ambient,
    (global_summary_wide$Global_Mean_ADM_Warmed - global_summary_wide$Global_Mean_ADM_Ambient) / global_summary_wide$Global_Mean_ADM_Ambient
  )
)

global_overall <- left_join(global_raw_overall, global_prop_overall, by = "Metric")
print(global_overall)

## 6. Plotting globally -------------------------------------------------------------------

# Global raw bar chart with SE bars
bar_col   <- "orange"
bar_width <- 0.5
err_width <- 0.15

  # Per-site raw deltas
  delta_by_site <- summary_table_all %>%
    tidyr::pivot_wider(
      id_cols = site,
      names_from = group,
      values_from = c(`Mean max`, `Mean ADM`)
    ) %>%
    dplyr::transmute(
      site,
      d_meanmax = `Mean max_Warmed` - `Mean max_Ambient`,
      d_adm     = `Mean ADM_Warmed`  - `Mean ADM_Ambient`
    )
  
  # Check; should be 11
  n_sites <- nrow(delta_by_site)
  
  # Plot
  global_bar_raw <- tibble::tibble(
    Metric     = c("Mean max", "Mean ADM"),
    mean_delta = c(mean(delta_by_site$d_meanmax, na.rm = TRUE),
                   mean(delta_by_site$d_adm,     na.rm = TRUE)),
    sd_delta   = c(sd(delta_by_site$d_meanmax,   na.rm = TRUE),
                   sd(delta_by_site$d_adm,       na.rm = TRUE))
  ) %>%
    dplyr::mutate(se_delta = sd_delta / sqrt(n_sites))
  
  p_global_raw_err <- ggplot(global_bar_raw, aes(x = Metric, y = mean_delta)) +
    geom_col(width = bar_width, fill = bar_col) +
    geom_errorbar(aes(ymin = mean_delta - se_delta,
                      ymax = mean_delta + se_delta),
                  width = err_width) +
    labs(title = "Raw temperature change",
         x = "Metric", y = "Temperature difference (°C)") +
    theme_minimal() +
    theme(
      # match pollution plots: faint vertical + horizontal majors
      panel.grid.major.x = element_line(color = "grey85"),
      panel.grid.major.y = element_line(color = "grey85"),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", hjust = 0.5),
      legend.position    = "none"
    )
  
  p_global_raw_err

# Global raw vs proportional differences heat map
  
  # Pivot summary_table_all so that each site has separate columns for Ambient and Warmed
  global_raw <- summary_table_all %>%
    pivot_wider(
      id_cols = site,
      names_from = group,
      values_from = c(`Mean max`, `Mean ADM`)
    ) %>%
    mutate(
      Diff_Mean_max = `Mean max_Warmed` - `Mean max_Ambient`,
      Diff_Mean_ADM  = `Mean ADM_Warmed`  - `Mean ADM_Ambient`
    )
  
  # Reshape raw differences into long format
  global_raw_long <- global_raw %>%
    dplyr::select(site, Diff_Mean_max, Diff_Mean_ADM) %>%
    tidyr::pivot_longer(
      cols      = c(Diff_Mean_max, Diff_Mean_ADM),
      names_to  = "Metric",
      values_to = "RawDiff"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        "Diff_Mean_max" = "Mean max",
        "Diff_Mean_ADM" = "Mean ADM"
      ),
      site = factor(site, levels = sort(unique(site)))
    )
  
  # Raw plot
  raw_heat <- ggplot(global_raw_long, aes(x = Metric, y = site, fill = RawDiff)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%+.1f°C", RawDiff)), 
              color = "black", size = 3) +
    scale_fill_gradient(
      low = "white", high = "red",
      name = "Raw increase (°C)",
      breaks = seq(
        0,
        ceiling(max(global_raw_long$RawDiff, na.rm = TRUE) * 2) / 2,
        by = 0.5
      ),
      labels = function(x) paste0(x, "°C")
    ) +
    labs(
      title = "Raw temperature change",
      x = "Metric", y = "Site"
    ) +
    scale_y_discrete(limits = rev(levels(global_raw_long$site))) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Global proportional differences
  
  # Pivot summary_table_all into wide format
  global_prop <- summary_table_all %>%
    pivot_wider(
      id_cols = site,
      names_from = group,
      values_from = c(`Mean max`, `Mean ADM`)
    ) %>%
    mutate(
      PropDiff_Mean_max = (`Mean max_Warmed` - `Mean max_Ambient`) / `Mean max_Ambient`,
      PropDiff_Mean_ADM = (`Mean ADM_Warmed` - `Mean ADM_Ambient`) / `Mean ADM_Ambient`
    )
  
  # Reshape proportional differences into long format
  global_prop_long <- global_prop %>%
    dplyr::select(site, PropDiff_Mean_max, PropDiff_Mean_ADM) %>%
    tidyr::pivot_longer(
      cols      = c(PropDiff_Mean_max, PropDiff_Mean_ADM),
      names_to  = "Metric",
      values_to = "PropDiff"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        "PropDiff_Mean_max" = "Mean max",
        "PropDiff_Mean_ADM"  = "Mean ADM"
      ),
      site = factor(site, levels = sort(unique(site)))
    )

  # Proportional plot
  prop_heat <- ggplot(global_prop_long, aes(x = Metric, y = site, fill = PropDiff)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%+.1f%%", PropDiff * 100)),
              color = "black", size = 3) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0,
      name = "Proportional\nincrease (%)",
      labels = scales::percent,
      breaks = seq(0, 0.10, by = 0.02)
    ) +
    labs(
      title = "Proportional temperature change",
      x = "Metric", y = "Site"
    ) +
    scale_y_discrete(limits = rev(levels(global_prop_long$site))) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
final_heat_panel <- raw_heat | prop_heat
final_heat_panel