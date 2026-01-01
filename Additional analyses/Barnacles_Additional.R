###############################################################################
# Barnacles (additional analyses)
#
#   Raw exploratory plots
#   Diagnostics + alternative model variants / transformations
#   Random-slope structure selection evidence (AIC + singularity flags)
#   Pseudo-count sensitivity (log transform evidence)
#   Baseline/intercept BLUPs + baseline deviation plots
#   Per-site stressor-effects “caterpillar” on the % scale
#
#   *This script may print warnings/errors that are informative for troubleshooting purposes.
###############################################################################

## 0. Setup -------------------------------------------------------------------
library(here)
library(tidyverse)
library(lubridate)
library(lmerTest)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom.mixed)
library(broom)
library(patchwork)
library(grid)
library(ggplot2)
library(scales)
library(ggtext)
library(ggnewscale)
library(DT)

set.seed(1)

## 0.1 Global constants

# Shared palette
TREAT_FILL <- c(
  "NP_Ambient" = "grey75",
  "NP_Warmed"  = "orange",
  "P_Ambient"  = "#006400",
  "P_Warmed"   = "#CC79A7"
)

## 0.2 Helper functions

# Harmonise site naming for joins/plots only
clean_site <- function(x) {
  dplyr::case_when(
    x == "Portugal1"        ~ "Portugal 1",
    x == "Portugal2"        ~ "Portugal 2",
    x == "Vietnam1"         ~ "Vietnam",
    x == "Vietnam"          ~ "Vietnam",
    x == "NZ"               ~ "New Zealand",
    x == "United Kingdom"   ~ "UK",
    TRUE                    ~ x
  )
}

# Clean labels used in caterpillar plots
pretty_site_labels <- function(x) {
  dplyr::recode(
    x,
    "Portugal 1" = "Portugal (Madeira)",
    "Portugal 2" = "Portugal (Porto)",
    "Vietnam"    = "Vietnam",
    .default     = x
  )
}

# Variance-covariance arrays can contain tiny negative values due to floating point rounding.
# We clamp <0 to 0 before sqrt to avoid NaNs.
safe_sqrt <- function(x) {
  v <- as.numeric(x)
  v[!is.finite(v) | v < 0] <- 0
  sqrt(v)
}

# Pearson-style diagnostics for lmer models:
#   residuals vs fitted
#   QQ plot
#   histogram + normal overlay
#   scale-location
make_diagnostics_lmer <- function(mod, resp_name) {
  mf <- model.frame(mod)
  mf$.fitted  <- fitted(mod)
  mf$.resid   <- resid(mod)
  mf$abs_sqrt <- sqrt(abs(mf$.resid))
  
  p1 <- ggplot(mf, aes(.fitted, .resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste0("Residuals vs fitted (", resp_name, ")"),
      x = "Fitted", y = "Residual"
    ) +
    theme_minimal()
  
  p2 <- ggplot(mf, aes(sample = .resid)) +
    stat_qq(alpha = 0.6) +
    stat_qq_line() +
    labs(title = paste0("QQ plot (", resp_name, ")")) +
    theme_minimal()
  
  p3 <- ggplot(mf, aes(.resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray80") +
    stat_function(
      fun  = dnorm,
      args = list(mean = mean(mf$.resid), sd = sd(mf$.resid))
    ) +
    labs(title = paste0("Residual histogram (", resp_name, ")")) +
    theme_minimal()
  
  p4 <- ggplot(mf, aes(.fitted, abs_sqrt)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste0("Scale–location (", resp_name, ")"),
      x = "Fitted", y = "√|Residual|"
    ) +
    theme_minimal()
  
  (p1 | p2) / (p3 | p4)
}

add_treatment_factor <- function(df) {
  df %>%
    mutate(
      Treatment = case_when(
        Pollution == 0 & Warming == 0 ~ "Non-Polluted + Ambient",
        Pollution == 0 & Warming == 1 ~ "Non-Polluted + Warm",
        Pollution == 1 & Warming == 0 ~ "Polluted + Ambient",
        Pollution == 1 & Warming == 1 ~ "Polluted + Warm"
      ) %>%
        factor(levels = c(
          "Non-Polluted + Ambient", "Non-Polluted + Warm",
          "Polluted + Ambient", "Polluted + Warm"
        ))
    )
}

## 1. Data preparation & aggregation ------------------------------------------

## 1.1 Load full dataset
data_raw <- readr::read_csv(
  here::here("Dataframes", "Full global HotMess data.csv"),
  show_col_types = FALSE
)

## 1.2 Recode treatments
data_prep <- data_raw %>%
  mutate(
    Warming   = if_else(`0-Pad colour` == "B", 1L, 0L),
    Pollution = if_else(`0-Pad region` == "P", 1L, 0L)
  )

## 1.3 Aggregate to pad-level summer means
pad_summary <- data_prep %>%
  group_by(`0-Country/site`, `0-Pad ID`) %>%
  summarise(
    `0-Pad region`   = first(`0-Pad region`),
    `0-Pad colour`   = first(`0-Pad colour`),
    `0-Latitude`     = first(`0-Latitude`),
    `4-Is_temperate` = first(`4-Is_temperate`),
    Warming          = first(Warming),
    Pollution        = first(Pollution),
    AvgBarnacles     = mean(`3-Barnacle total count`, na.rm = TRUE),
    TotalBarnacles   = sum(`3-Barnacle total count`,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(AvgBarnacles))

## 1.4 Drop 'zero-barnacle' sites (mean == 0 across summer)
zero_barnacle_sites <- pad_summary %>%
  group_by(`0-Country/site`) %>%
  summarise(mean_barn_site = mean(AvgBarnacles, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_barn_site == 0) %>%
  pull(`0-Country/site`)

## 1.5 Build model response + plot keys
pad_summary <- pad_summary %>%
  filter(!(`0-Country/site` %in% zero_barnacle_sites)) %>%
  mutate(
    # Transform candidates used for evidence in Section 2
    log_AvgBarnacles   = log(AvgBarnacles + 1),
    sqrt_AvgBarnacles  = sqrt(AvgBarnacles + 0.5),
    asinh_AvgBarnacles = asinh(AvgBarnacles),
    
    # Clean site label for plotting/joins only
    Site_clean = clean_site(`0-Country/site`),
    
    # Group label for quick residual checks
    group  = interaction(Warming, Pollution, sep = "×"),
  )

## 1.6 Exploratory raw visualisation
pad_summary <- add_treatment_factor(pad_summary)

TREAT_LABEL_COLS <- c(
  "Non-Polluted + Ambient" = unname(TREAT_FILL["NP_Ambient"]),
  "Non-Polluted + Warm"    = unname(TREAT_FILL["NP_Warmed"]),
  "Polluted + Ambient"     = unname(TREAT_FILL["P_Ambient"]),
  "Polluted + Warm"        = unname(TREAT_FILL["P_Warmed"])
)

raw_plot <- ggplot(pad_summary, aes(x = Treatment, y = AvgBarnacles, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(fill = Treatment),
    shape  = 21,
    colour = "black",
    stroke = 0.5,
    width  = 0.15,
    size   = 2,
    alpha  = 0.6
  ) +
  scale_fill_manual(values = TREAT_LABEL_COLS) +
  labs(
    title = "Barnacle abundance by Warming × Pollution treatment",
    x = NULL,
    y = "Mean barnacle count (pad-level summer mean)"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(raw_plot)

## 2. Model variants & diagnostics --------------------------------------------
#   Summarise evidence for the modelling decisions used in main code.

## 2.1 Baseline: raw response + full random slopes (evidence)
barn_fullRE_raw <- lmerTest::lmer(
  AvgBarnacles ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  data = pad_summary,
  REML = FALSE
) # singular

print(summary(barn_fullRE_raw))
print(make_diagnostics_lmer(barn_fullRE_raw, "Gaussian LMM (AvgBarnacles)"))

## 2.2 Transformations: sqrt and log (full random slopes)
barn_fullRE_sqrt <- update(barn_fullRE_raw, sqrt_AvgBarnacles ~ .) # singular
barn_fullRE_log  <- update(barn_fullRE_raw, log_AvgBarnacles  ~ .)

print(summary(barn_fullRE_sqrt))
print(summary(barn_fullRE_log))

print(make_diagnostics_lmer(barn_fullRE_sqrt, "sqrt(Avg+0.5)"))
print(make_diagnostics_lmer(barn_fullRE_log,  "log(Avg+1)"))

## 2.3 Gamma alternatives (glmmTMB)
# Barnacle counts are not truly 'continuous', so Gamma is evidence only here as a troubleshooting metric.
mod_gamma_id <- glmmTMB(
  AvgBarnacles + 0.1 ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  family = Gamma(link = "identity"),
  data   = pad_summary
) # several warnings

mod_gamma_log <- glmmTMB(
  AvgBarnacles + 0.1 ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  family = Gamma(link = "log"),
  data   = pad_summary
) # convergence issues

print(summary(mod_gamma_id))
print(summary(mod_gamma_log))

## 2.4 Random-slope structure selection (AIC + singularity evidence)
slope_sets <- list(
  intercept_only     = character(0),
  warming            = "Warming",
  pollution          = "Pollution",
  interaction        = "Warming:Pollution",
  W_plus_P           = c("Warming", "Pollution"),
  W_plus_INT         = c("Warming", "Warming:Pollution"),
  P_plus_INT         = c("Pollution", "Warming:Pollution"),
  full_random_slopes = c("Warming", "Pollution", "Warming:Pollution")
)

make_rand <- function(slopes) {
  if (length(slopes) == 0) return("(1 | `0-Country/site`)")
  paste0("(1 + ", paste(slopes, collapse = " + "), " | `0-Country/site`)")
}

re_selection_tbl <- purrr::imap_dfr(slope_sets, ~{
  rand_chunk <- make_rand(.x)
  fmla <- as.formula(paste0("log_AvgBarnacles ~ Warming * Pollution + ", rand_chunk))
  
  mod <- try(lmerTest::lmer(fmla, data = pad_summary, REML = FALSE), silent = TRUE)
  if (inherits(mod, "try-error")) {
    return(tibble(config = .y, singular = NA, AIC = NA))
  }
  
  tibble(
    config   = .y,
    singular = lme4::isSingular(mod, tol = 1e-4),
    AIC      = AIC(mod)
  )
})

re_selection_tbl %>%
  arrange(singular, AIC) %>%
  mutate(allowed = if_else(singular == FALSE, "Y", "N")) %>%
  select(config, allowed, AIC) %>%
  print() # Full slopes supported, and relatively little difference in AIC comparatively across 'best' models

## 2.5 Pseudo-count sensitivity for log transform model
# Main code uses log(Avg+1). Here we demonstrate qualitative consistency with variation in pseudocounts.
pseudos <- c(0.0625, 0.1, 0.125, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

sens_all <- purrr::map_dfr(pseudos, function(k) {
  df2 <- pad_summary %>% mutate(logB = log(AvgBarnacles + k))
  
  mod_k <- lmerTest::lmer(
    logB ~ Warming * Pollution +
      (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
    data = df2,
    REML = FALSE
  )
  
  fe_k <- broom.mixed::tidy(mod_k, effects = "fixed") %>%
    filter(term %in% c("Warming", "Pollution", "Warming:Pollution")) %>%
    select(term, estimate, p.value) %>%
    pivot_wider(
      names_from = term,
      values_from = c(estimate, p.value),
      names_sep = "_"
    )
  
  tibble(pseudo = k, AIC = AIC(mod_k)) %>% bind_cols(fe_k)
})

print(sens_all) # effect sizes and significance relatively stable across any choice; retain +1 for interpretability.

## 3. Site-level insights: intercept BLUPs --------------------------
model_final <- lmerTest::lmer(
  log_AvgBarnacles ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  data = pad_summary,
  REML = FALSE
)

re_block <- ranef(model_final, condVar = TRUE)$`0-Country/site`
postVar  <- attr(re_block, "postVar")
re_terms <- colnames(re_block)

# Helper to pull var-cov entry for a given term pair (i, j) across sites
pv_ij <- function(term_i, term_j) {
  ii <- which(re_terms == term_i)
  jj <- which(re_terms == term_j)
  postVar[ii, jj, ]
}

site_blups <- tibble(
  Site     = rownames(re_block),
  BLUP_log = re_block[, "(Intercept)"],
  se_log   = safe_sqrt(pv_ij("(Intercept)", "(Intercept)"))
) %>%
  mutate(Site = clean_site(Site))

alpha_log    <- unname(fixef(model_final)["(Intercept)"])
global_count <- exp(alpha_log)

# Relative baseline deviation (%)
site_pct <- site_blups %>%
  mutate(
    ratio     = exp(BLUP_log),
    pct_dev   = (ratio - 1) * 100,
    lower_pct = (exp(BLUP_log - 1.96 * se_log) - 1) * 100,
    upper_pct = (exp(BLUP_log + 1.96 * se_log) - 1) * 100
  )

# Absolute baseline deviation (counts)
site_raw <- site_blups %>%
  mutate(
    site_count     = exp(alpha_log + BLUP_log),
    lower_count    = exp(alpha_log + BLUP_log - 1.96 * se_log),
    upper_count    = exp(alpha_log + BLUP_log + 1.96 * se_log),
    raw_dev        = site_count - global_count,
    lower_raw_dev  = lower_count - global_count,
    upper_raw_dev  = upper_count - global_count
  )

p_pct <- ggplot(site_pct, aes(x = pct_dev, y = reorder(Site, pct_dev))) +
  geom_errorbarh(aes(xmin = lower_pct, xmax = upper_pct),
                 height = 0.2, color = "gray60") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Baseline % deviation from global mean",
       x = "% deviation", y = NULL) +
  theme_minimal()

p_raw <- ggplot(site_raw, aes(x = raw_dev, y = reorder(Site, raw_dev))) +
  geom_errorbarh(aes(xmin = lower_raw_dev, xmax = upper_raw_dev),
                 height = 0.2, color = "gray60") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Baseline raw count deviation from global mean",
       x = "Count deviation", y = NULL) +
  theme_minimal()

baseline_cat_plot <- p_pct + p_raw +
  plot_annotation(
    title    = "Site-level baseline barnacle abundance deviations",
    subtitle = "Left = relative % deviation; Right = absolute deviation in mean count"
  )

print(baseline_cat_plot)

## 4. Site-level stressor effects (caterpillar plot on % scale; main code shows on model scale) --------------------
re_site  <- ranef(model_final, condVar = TRUE)$`0-Country/site`
pv_array <- attr(re_site, "postVar")
re_terms <- colnames(re_site)
fix_eff  <- fixef(model_final)

pv_ij2 <- function(term_i, term_j) {
  ii <- which(re_terms == term_i)
  jj <- which(re_terms == term_j)
  pv_array[ii, jj, ]
}

# Site ordering north to south
site_order <- pad_summary %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site_raw = `0-Country/site`,
    Site     = clean_site(`0-Country/site`),
    Latitude = as.numeric(`0-Latitude`)
  ) %>%
  arrange(desc(Latitude)) %>%
  pull(Site)

y_levels <- rev(c("<b>Global fixed effect</b>", site_order))

site_effects <- tibble(
  Site_raw = rownames(re_site),
  Site     = clean_site(rownames(re_site)),
  u_W      = re_site[, "Warming"],
  u_P      = re_site[, "Pollution"],
  u_INT    = re_site[, "Warming:Pollution"]
) %>%
  mutate(
    warm_log = unname(fix_eff["Warming"]) + u_W,
    poll_log = unname(fix_eff["Pollution"]) + u_P,
    int_log  = unname(fix_eff["Warming:Pollution"]) + u_INT,
    
    warm_se  = safe_sqrt(pv_ij2("Warming", "Warming")),
    poll_se  = safe_sqrt(pv_ij2("Pollution", "Pollution")),
    int_se   = safe_sqrt(pv_ij2("Warming:Pollution", "Warming:Pollution"))
  )

all_df <- bind_rows(
  site_effects %>% transmute(Site, Term = "Warming",     slope_log = warm_log, se_log = warm_se),
  site_effects %>% transmute(Site, Term = "Pollution",   slope_log = poll_log, se_log = poll_se),
  site_effects %>% transmute(Site, Term = "Interaction", slope_log = int_log,  se_log = int_se)
) %>%
  mutate(
    lower_log = slope_log - 1.96 * se_log,
    upper_log = slope_log + 1.96 * se_log,
    pct       = 100 * (exp(slope_log) - 1),
    lower_pct = 100 * (exp(lower_log) - 1),
    upper_pct = 100 * (exp(upper_log) - 1),
    Term      = factor(Term, levels = c("Warming", "Pollution", "Interaction")),
    Site      = factor(as.character(Site), levels = y_levels)
  )

vc_fix <- as.matrix(vcov(model_final))

global_df <- tibble(
  Term = factor(c("Warming", "Pollution", "Interaction"),
                levels = c("Warming", "Pollution", "Interaction")),
  slope_log = c(unname(fix_eff["Warming"]),
                unname(fix_eff["Pollution"]),
                unname(fix_eff["Warming:Pollution"])),
  se_log = c(
    safe_sqrt(vc_fix["Warming", "Warming"]),
    safe_sqrt(vc_fix["Pollution", "Pollution"]),
    safe_sqrt(vc_fix["Warming:Pollution", "Warming:Pollution"])
  )
) %>%
  mutate(
    lower_log = slope_log - 1.96 * se_log,
    upper_log = slope_log + 1.96 * se_log,
    pct       = 100 * (exp(slope_log) - 1),
    lower_pct = 100 * (exp(lower_log) - 1),
    upper_pct = 100 * (exp(upper_log) - 1),
    Site      = factor("<b>Global fixed effect</b>", levels = y_levels)
  )

p_pct_barn <- ggplot() +
  geom_rect(
    data = global_df,
    aes(xmin = lower_pct, xmax = upper_pct),
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.3
  ) +
  geom_vline(
    data = global_df,
    aes(xintercept = pct),
    linetype = "dotted",
    colour = "black"
  ) +
  geom_point(
    data = global_df,
    aes(x = pct, y = Site, fill = Term),
    shape = 23, size = 5, colour = "black"
  ) +
  geom_errorbarh(
    data = all_df,
    aes(y = Site, xmin = lower_pct, xmax = upper_pct),
    height = 0.2, colour = "#666666"
  ) +
  geom_point(
    data = all_df,
    aes(x = pct, y = Site, colour = Term),
    size = 3
  ) +
  scale_color_manual(
    values = c(Warming = "orange", Pollution = "darkgreen", Interaction = "#CC79A7"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(Warming = "orange", Pollution = "darkgreen", Interaction = "#CC79A7"),
    guide = "none"
  ) +
  facet_wrap(
    ~Term,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(
      Warming = "Warming effect",
      Pollution = "Pollution effect",
      Interaction = "Interaction effect"
    )
  ) +
  scale_x_continuous(name = "% change in mean barnacle count") +
  scale_y_discrete(
    limits = y_levels,
    labels = pretty_site_labels,
    name = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold"),
    panel.spacing.x = unit(1.5, "lines"),
    axis.text.y     = ggtext::element_markdown(size = 9),
    legend.position = "none"
  )

print(p_pct_barn)