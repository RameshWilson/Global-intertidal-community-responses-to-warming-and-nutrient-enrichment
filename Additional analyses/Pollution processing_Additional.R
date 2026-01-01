## 0. Setup -------------------------------------------------------------------

# install.packages("here")   # Install once per machine should code and data be read from GitHub
library(here)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(readr)
library(stringr)
library(readr)  
library(forcats)

input_csv <- here::here("Dataframes", "Global pollution raw data.csv")

# Numeric parser (tolerates dashes, 'nd', blanks, etc.)
to_num <- function(x) parse_number(
  as.character(x),
  na = c("NA", "NaN", "", "—", "–", "-", "nd", "ND")
)

raw <- read_csv(input_csv)

# Convert to long format
#   Site/region/replicate
#   Pollutant (Nitrate/Phosphate)
#   Value (mg/L)

long <- raw %>%
  rename(
    site           = Site,
    treatment      = Region,
    replicate      = Replicate,
    nitrate_mgL    = `Nitrate (mg/L)`,
    phosphate_mgL  = `Phosphate (mg/L)`
  ) %>%
  mutate(
    nitrate_mgL    = to_num(nitrate_mgL),
    phosphate_mgL  = to_num(phosphate_mgL)
  ) %>%
  pivot_longer(
    cols = c(nitrate_mgL, phosphate_mgL),
    names_to  = "pollutant",
    values_to = "value_mgL"
  ) %>%
  mutate(
    pollutant = dplyr::recode(
      as.character(pollutant),
      nitrate_mgL   = "Nitrate",
      phosphate_mgL = "Phosphate",
      .default = pollutant
    )
  ) %>%
  filter(!is.na(value_mgL))

# Check: should show 847 rows, 11 sites, and both pollutants
cat("Rows:", nrow(long), " | Sites:", dplyr::n_distinct(long$site),
    " | Pollutants:", paste(unique(long$pollutant), collapse = ", "), "\n")

long %>%
  count(site, treatment, pollutant, name = "n_rows") %>%
  arrange(pollutant, site, treatment) %>%
  print(n = Inf)

# Normalise if needed
long_clean <- long %>%
  mutate(
    site      = str_squish(site),
    treatment = str_squish(treatment),
    treatment = str_replace_all(treatment, "[\u2013\u2014]", "-"),  # en/em dash -> hyphen
    treatment = dplyr::recode(
      treatment,
      "Non polluted" = "Non-polluted",
      "non-polluted" = "Non-polluted",
      "non polluted" = "Non-polluted",
      "polluted"     = "Polluted",
      .default = treatment
    )
  )

# Site-level means by treatment
site_means <- long_clean %>%
  group_by(site, pollutant, treatment) %>%
  summarise(mean_mgL = mean(value_mgL, na.rm = TRUE), .groups = "drop")

# Pivot to wide (ambient vs polluted)
wide <- site_means %>%
  dplyr::select(site, pollutant, treatment, mean_mgL) %>%
  tidyr::pivot_wider(
    id_cols    = c(site, pollutant),
    names_from = treatment,
    values_from= mean_mgL
  ) %>%
  dplyr::rename(
    ambient_mean_mgL  = `Non-polluted`,
    polluted_mean_mgL = `Polluted`
  )

# Check
wide %>%
  dplyr::arrange(pollutant, site) %>%
  dplyr::select(site, pollutant, ambient_mean_mgL, polluted_mean_mgL) %>%
  dplyr::mutate(across(c(ambient_mean_mgL, polluted_mean_mgL), ~round(.x, 4))) %>%
  print(n = Inf)

missing_pairs <- wide %>%
  filter(is.na(ambient_mean_mgL) | is.na(polluted_mean_mgL)) %>%
  arrange(pollutant, site)

if (nrow(missing_pairs) > 0) {
  message("Heads-up: some site×pollutant are missing one side:")
  print(missing_pairs)
} else {
  message("All site×pollutant combos have both ambient and polluted means ✅")
}

## 1. Pseudocount and proportional increases -------------------------------------------------------------------

#   - First collapse each site × treatment × pollutant to a single summer mean (mg/L)
#   - To avoid divide-by-zero amplificatioins at oligotrophic sites, we apply a small, data-adaptive floor ('epsilon') to the denominator only when the ambient mean is incredibly small
#   - We then compute the site-level proportional increase as [polluted - ambient] divided by the larger of ambient or epsilon, and flag any site where the floor was actually used
#   - For global summaries, primary mean excludes any site where the floor was used (i.e., flagged sites with near-zero ambient), as those few cases can artificially dominate the average

# Controls (mg/L)
eps_min    <- 1e-4   # Absolute floor (negligible at mg/L scale)
eps_factor <- 0.05   # 5% of the 5th percentile of non-zero ambient means

# Data-driven per-pollutant (N/P) epsilon from lower tail of ambient means
eps_tbl <- wide %>%
  dplyr::group_by(pollutant) %>%
  dplyr::summarise(
    p5_nonzero = {
      v <- ambient_mean_mgL[is.finite(ambient_mean_mgL) & ambient_mean_mgL > 0]
      if (length(v) == 0) NA_real_ else as.numeric(stats::quantile(v, 0.05, na.rm = TRUE))
    },
    pseudocount_mgL = pmax(eps_min, eps_factor * tidyr::replace_na(p5_nonzero, 0)),
    .groups = "drop"
  )
print(eps_tbl)

# Proportional increases with stabilised denominator where needed
results <- wide %>%
  dplyr::left_join(eps_tbl, by = "pollutant") %>%
  dplyr::mutate(
    denom = dplyr::if_else(is.na(ambient_mean_mgL), NA_real_,
                           pmax(ambient_mean_mgL, pseudocount_mgL)),
    prop_increase     = (polluted_mean_mgL - ambient_mean_mgL) / denom,
    prop_increase_pct = 100 * prop_increase,
    denom_was_epsilon = !is.na(ambient_mean_mgL) & ambient_mean_mgL <= pseudocount_mgL
  ) %>%
  dplyr::arrange(pollutant, site)

# Check
results %>%
  dplyr::select(site, pollutant, ambient_mean_mgL, polluted_mean_mgL,
                pseudocount_mgL, prop_increase, denom_was_epsilon) %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4))) %>%
  print(n = Inf)
  # Only Portugal 2 Nitrate hit the epsilon floor (denom_was_epsilon = TRUE) because its ambient mean is 0 mg/L
  # Yields a huge proportional increase (≈ 5082×), which will inflate any plain mean.

## 2. Global summaries -------------------------------------------------------------------
#
# This computes cross-site stats two ways: 
#   - Inclusive of all sites
#   - Excluding sites where the epsilon floor was used (i.e., ambient ≈ 0) so a single site can’t dominate the mean

# Inclusive (all sites)
global_inclusive <- results %>%
  group_by(pollutant) %>%
  summarise(
    n_sites     = sum(!is.na(prop_increase)),
    mean_prop   = mean(prop_increase, na.rm = TRUE),
    sd_prop     = ifelse(n_sites > 1, sd(prop_increase, na.rm = TRUE), NA_real_),
    median_prop = median(prop_increase, na.rm = TRUE),
    q1_prop     = quantile(prop_increase, 0.25, na.rm = TRUE),
    q3_prop     = quantile(prop_increase, 0.75, na.rm = TRUE),
    mean_pct    = 100 * mean_prop,
    sd_pct      = 100 * sd_prop,
    .groups     = "drop"
  )

# Exclude sites where epsilon replaced the denominator (ambient ≈ 0)
global_excl_eps <- results %>%
  filter(!denom_was_epsilon) %>%
  group_by(pollutant) %>%
  summarise(
    n_sites_used     = sum(!is.na(prop_increase)),
    mean_prop_excl   = mean(prop_increase, na.rm = TRUE),
    sd_prop_excl     = ifelse(n_sites_used > 1, sd(prop_increase, na.rm = TRUE), NA_real_),
    median_prop_excl = median(prop_increase, na.rm = TRUE),
    q1_prop_excl     = quantile(prop_increase, 0.25, na.rm = TRUE),
    q3_prop_excl     = quantile(prop_increase, 0.75, na.rm = TRUE),
    mean_pct_excl    = 100 * mean_prop_excl,
    sd_pct_excl      = 100 * sd_prop_excl,
    .groups          = "drop"
  )

global_inclusive %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

global_excl_eps %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print(n = Inf)

## 3. Plotting -------------------------------------------------------------------

bar_col   <- "darkgreen"
bar_width <- 0.5
err_width <- 0.15

logger_bar_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_line(color = "grey85"),
      panel.grid.major.y = element_line(color = "grey85"),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", hjust = 0.5),
      legend.position    = "none"
    )
}

pct_scale <- scale_y_continuous(
  labels = label_number(accuracy = 1, big.mark = ",", suffix = "%"),
  expand = expansion(mult = c(0.02, 0.15))   # headroom for labels if present
)

# Ensure consistent order
lvl_pol <- c("Nitrate","Phosphate")

# Mean ± SE
plot_global_mean_se <- global_excl_eps %>%
  mutate(
    pollutant      = factor(pollutant, levels = lvl_pol),
    mean_pct_excl  = 100 * mean_prop_excl,
    sd_pct_excl    = 100 * sd_prop_excl,
    se_pct_excl    = sd_pct_excl / sqrt(n_sites_used)
  ) %>%
  ggplot(aes(x = pollutant, y = mean_pct_excl)) +
  geom_col(width = bar_width, fill = bar_col) +
  geom_errorbar(aes(ymin = mean_pct_excl - se_pct_excl,
                    ymax = mean_pct_excl + se_pct_excl),
                width = err_width) +
  pct_scale +
  labs(
    title = "Global proportional increase",
    subtitle = "Mean ± SE across sites (excluding epsilon-denominator sites)",
    x = "Pollutant", y = "Proportional increase (%)"
  ) +
  logger_bar_theme()

print(plot_global_mean_se)

# Median + IQR
plot_global_median_iqr <- global_excl_eps %>%
  mutate(
    pollutant = factor(pollutant, levels = lvl_pol),
    med_pct   = 100 * median_prop_excl,
    q1_pct    = 100 * q1_prop_excl,
    q3_pct    = 100 * q3_prop_excl
  ) %>%
  ggplot(aes(x = pollutant, y = med_pct)) +
  geom_col(width = bar_width, fill = bar_col) +
  geom_errorbar(aes(ymin = q1_pct, ymax = q3_pct), width = err_width) +
  pct_scale +
  labs(
    title = "Global proportional increase",
    subtitle = "Median with IQR across sites (excluding epsilon-denominator sites)",
    x = "Pollutant", y = "Proportional increase (%)"
  ) +
  logger_bar_theme()

print(plot_global_median_iqr)

# Mean ± SE + value labels
plot_global_mean_se_lab <- global_excl_eps %>%
  mutate(
    pollutant      = factor(pollutant, levels = lvl_pol),
    mean_pct_excl  = 100 * mean_prop_excl,
    sd_pct_excl    = 100 * sd_prop_excl,
    se_pct_excl    = sd_pct_excl / sqrt(n_sites_used),
    label_mean     = number(mean_pct_excl, accuracy = 1, big.mark = ",", suffix = "%")
  ) %>%
  ggplot(aes(x = pollutant, y = mean_pct_excl)) +
  geom_col(width = bar_width, fill = bar_col) +
  geom_errorbar(aes(ymin = mean_pct_excl - se_pct_excl,
                    ymax = mean_pct_excl + se_pct_excl),
                width = err_width) +
  geom_text(aes(label = label_mean), vjust = -0.5, size = 4) +
  pct_scale +
  labs(
    title = "Global proportional increase",
    subtitle = "Mean ± SE across sites (excluding epsilon-denominator sites)",
    x = "Pollutant", y = "Proportional increase (%)"
  ) +
  logger_bar_theme() +
  coord_cartesian(clip = "off")   # keep labels visible

print(plot_global_mean_se_lab)

# Median + IQR + value labels
plot_global_median_iqr_lab <- global_excl_eps %>%
  mutate(
    pollutant   = factor(pollutant, levels = lvl_pol),
    med_pct     = 100 * median_prop_excl,
    q1_pct      = 100 * q1_prop_excl,
    q3_pct      = 100 * q3_prop_excl,
    label_median= number(med_pct, accuracy = 1, big.mark = ",", suffix = "%")
  ) %>%
  ggplot(aes(x = pollutant, y = med_pct)) +
  geom_col(width = bar_width, fill = bar_col) +
  geom_errorbar(aes(ymin = q1_pct, ymax = q3_pct), width = err_width) +
  geom_text(aes(label = label_median), vjust = -0.5, size = 4) +
  pct_scale +
  labs(
    title = "Global proportional increase",
    subtitle = "Median with IQR across sites (excluding epsilon-denominator sites)",
    x = "Pollutant", y = "Proportional increase (%)"
  ) +
  logger_bar_theme() +
  coord_cartesian(clip = "off")

print(plot_global_median_iqr_lab)

# Heat map
bar_col <- "darkgreen"
p_cap   <- 0.90  # <-- use 0.90 for 90th percentile cap

heat_df <- results %>%
  mutate(
    pollutant = factor(pollutant, levels = c("Nitrate","Phosphate")),
    site = fct_reorder(site, prop_increase_pct, .fun = mean, .na_rm = TRUE)
  )

cap <- quantile(heat_df$prop_increase_pct[!heat_df$denom_was_epsilon], p_cap, na.rm = TRUE)

heat_df <- heat_df %>%
  mutate(
    fill_val  = pmin(prop_increase_pct, cap),
    label_txt = ifelse(is.na(prop_increase_pct), "",
                       paste0(number(prop_increase_pct, accuracy = 0.1),
                              ifelse(denom_was_epsilon, " ★", ""))),
    label_col = ifelse(fill_val >= 0.85 * cap, "white", "black")  # tweak contrast threshold a touch
  )

p_heat <- ggplot(heat_df, aes(x = pollutant, y = site, fill = fill_val)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = bar_col,
                      limits = c(0, cap), oob = squish,
                      name = "Prop. increase (%)",
                      labels = label_number(accuracy = 0.1)) +
  geom_text(aes(label = label_txt, color = label_col), size = 3) +
  scale_color_identity() +
  labs(title = "Per-site proportional increase",
       subtitle = paste0("Single-hue scale capped at the ",
                         percent(p_cap, accuracy = 1), " of non-ε sites (cap = ",
                         number(cap, accuracy = 1), "%); ★ = epsilon used"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

plot(p_heat)