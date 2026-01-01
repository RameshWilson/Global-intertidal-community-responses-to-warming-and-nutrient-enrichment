###############################################################################
# Barnacle abundance (Main analyses)
#
#   Plate-level summer means (AvgBarnacles)
#   Final Gaussian LMM on log(AvgBarnacles + 1)
#   * Random effects: (1 + Warming + Pollution + Warming:Pollution | Site)
#   Fixed effects + EMM contrasts
#   Interaction classification
#   Conditional plot
#   Site-level stressor responses (log-scale caterpillar BLUP plot + biogeography tables)
#   Meta-regressions (warming sensitivity; pollution sensitivity)
#   Diagnostics & influence key pollution meta-regression
###############################################################################

## 0. Setup -------------------------------------------------------------------
library(here)
library(tidyverse)
library(lmerTest)
library(lme4)
library(emmeans)
library(broom.mixed)
library(broom)
library(patchwork)
library(ggtext)
library(grid)
library(scales)

set.seed(1)

## 0.1 Global constants --------------------------------------------------------

# Treatment fills
TREAT_FILL <- c(
  "NP_Ambient" = "grey75",
  "NP_Warmed"  = "orange",
  "P_Ambient"  = "#006400",
  "P_Warmed"   = "#CC79A7"
)

# Meta-regression strata palettes
PAL_STRATA <- c("Non-Polluted" = "grey75", "Polluted" = "#006400")
PAL_WARM   <- c("Ambient" = "grey75", "Warmed" = "orange")

# Interaction equivalence band on exponentiated interaction term (±5%)
# For log(Avg+1) Gaussian LMM, exp(beta_int) is an interaction ratio.
EQUIV_LO <- 0.95
EQUIV_HI <- 1.05

# Pseudocount used for log10 plotting
PSEUDO_COUNT <- 0.1

## 0.2 Helper functions --------------------------------------------------------

# emmeans/broom outputs use several different CI column conventions
#   (lower.CL / upper.CL, conf.low / conf.high, etc.).
# This helper standardises to `lower` and `upper` so downstream code is stable.
standardize_ci <- function(df) {
  lower_candidates <- c("lower.CL", "asymp.LCL", "LCL", "conf.low")
  upper_candidates <- c("upper.CL", "asymp.UCL", "UCL", "conf.high")
  
  lower_name <- lower_candidates[lower_candidates %in% names(df)][1]
  upper_name <- upper_candidates[upper_candidates %in% names(df)][1]
  
  if (is.na(lower_name) || is.na(upper_name)) {
    stop("standardize_ci(): Could not find CI columns. Present: ",
         paste(names(df), collapse = ", "))
  }
  
  df %>% dplyr::rename(lower = !!lower_name, upper = !!upper_name)
}

# lme4 posterior variance arrays can contain tiny negative values due to rounding.
# We clamp <0 to 0 before sqrt to avoid NaNs.
safe_sqrt <- function(x) {
  v <- as.numeric(x)
  v[!is.finite(v) | v < 0] <- 0
  sqrt(v)
}

# Harmonise site naming for joins/plots only.
clean_site <- function(x) {
  dplyr::case_when(
    x == "Portugal1" ~ "Portugal 1",
    x == "Portugal2" ~ "Portugal 2",
    x == "Vietnam1"  ~ "Vietnam",
    x == "Vietnam"   ~ "Vietnam",
    x == "NZ"        ~ "New Zealand",
    x == "United Kingdom" ~ "UK",
    TRUE ~ x
  )
}

# Clean labels used in caterpillar plots (names for display only)
pretty_site_labels <- function(x) {
  dplyr::recode(
    x,
    "Portugal 1" = "Portugal (Madeira)",
    "Portugal 2" = "Portugal (Porto)",
    "Vietnam"    = "Vietnam",
    .default     = x
  )
}

## 1. Data preparation & aggregation ------------------------------------------

## 1.1 Load full dataset
data_raw <- readr::read_csv(
  here::here("Dataframes", "Full global HotMess data.csv"),
  show_col_types = FALSE
)

## 1.2 Recode treatments + build total grazer count
# Warming:   Black pad = warmed (1), White pad = ambient (0)
# Pollution: P region  = polluted (1), N region = non-polluted (0)
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
  filter(is.finite(AvgBarnacles))

## 1.4 Drop 'zero-barnacle' sites (mean == 0 across summer)
# These sites contribute no information about proportional change in barnacle abundance.
zero_barnacle_sites <- pad_summary %>%
  group_by(`0-Country/site`) %>%
  summarise(mean_barn_site = mean(AvgBarnacles, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_barn_site == 0) %>%
  pull(`0-Country/site`)
print(zero_barnacle_sites) # Argentina

pad_summary <- pad_summary %>%
  filter(!(`0-Country/site` %in% zero_barnacle_sites))

## 1.5 Build model response + plot keys
pad_summary <- pad_summary %>%
  mutate(
    log_AvgBarnacles = log(AvgBarnacles + 1),
    Site_clean       = clean_site(`0-Country/site`)
  )

## 2. Final model specification -----------------------------------------------

## 2.1 Final LMM on log(AvgBarnacles + 1) with full random-effects structure
model_final <- lmerTest::lmer(
  log_AvgBarnacles ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  data = pad_summary,
  REML = FALSE
)

summary(model_final)

## 2.2 Quick model checks (numeric only)
# Fuller diagnostics / structure-selection evidence live in ADDITIONAL.
conv_msgs <- tryCatch(model_final@optinfo$conv$lme4$messages, error = function(e) NULL)

model_check_tbl <- tibble(
  singular      = lme4::isSingular(model_final, tol = 1e-4),
  n_conv_msgs   = if (is.null(conv_msgs)) 0L else length(conv_msgs)
)

print(model_check_tbl)
if (!is.null(conv_msgs)) print(conv_msgs)
# Supported random effects structure and clean convergence

## 3. Main results ------------------------------------------------------------

## 3.1 Fixed effects table (log scale + ratios + % change)
fe_tbl <- broom.mixed::tidy(
  model_final,
  effects  = "fixed",
  conf.int = TRUE
) %>%
  rename(df_obs = df) %>%
  transmute(
    term,
    estimate_log = estimate,
    std_error    = std.error,
    t_value      = statistic,
    df           = df_obs,
    p_value      = p.value,
    ratio        = exp(estimate),
    ratio_LCL    = exp(conf.low),
    ratio_UCL    = exp(conf.high),
    pct_change   = 100 * (ratio - 1)
  )

print(fe_tbl)

## 3.2 EMM contrasts
emm_link <- emmeans(model_final, ~ Warming * Pollution, type = "link")

# Warmed vs Ambient within each Pollution stratum
warm_contr <- contrast(emm_link, method = "revpairwise", by = "Pollution") %>%
  summary(infer = TRUE) %>%
  standardize_ci() %>%
  transmute(
    Contrast = "Warmed vs Ambient",
    Stratum  = if_else(Pollution == 0, "Non-Polluted", "Polluted"),
    Ratio    = exp(estimate),
    LCL      = exp(lower),
    UCL      = exp(upper),
    `Δ %`    = 100 * (Ratio - 1),
    p        = p.value
  )

# Polluted vs Non-polluted within each Warming stratum
poll_contr <- contrast(emm_link, method = "revpairwise", by = "Warming") %>%
  summary(infer = TRUE) %>%
  standardize_ci() %>%
  transmute(
    Contrast = "Polluted vs Non-polluted",
    Stratum  = if_else(Warming == 0, "Ambient", "Warmed"),
    Ratio    = exp(estimate),
    LCL      = exp(lower),
    UCL      = exp(upper),
    `Δ %`    = 100 * (Ratio - 1),
    p        = p.value
  )

emm_table_barnacles <- bind_rows(warm_contr, poll_contr)
print(emm_table_barnacles)

## 3.3 Interaction classification (±5% ratio band + p-value gate)

# We classify the Warming×Pollution interaction on the model’s natural (log) scale,
# but report exp(beta_int) as an interaction ratio using:
#   - p<0.05 gate on the interaction term
#   - ±5% equivalence band on ratio: [0.95, 1.05]
# Directional labels use additive vs observed changes on the log scale:
#   d_add = beta_W + beta_P
#   d_obs = d_add + beta_WP
extract_barnacle_interaction_class <- function(model,
                                               alpha    = 0.05,
                                               equiv_lo = EQUIV_LO,
                                               equiv_hi = EQUIV_HI,
                                               label    = "Whole Season") {
  
  fe <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  # If CI columns are missing for any reason, fall back to Wald
  if (!all(c("conf.low", "conf.high") %in% names(fe))) {
    fe <- fe %>%
      mutate(
        conf.low  = estimate - 1.96 * std.error,
        conf.high = estimate + 1.96 * std.error
      )
  }
  
  get_val <- function(term_name, col_name) {
    out <- fe %>%
      filter(.data$term == term_name) %>%
      slice(1) %>%
      pull(.data[[col_name]])
    if (length(out) == 0) NA_real_ else out
  }
  
  b0  <- get_val("(Intercept)", "estimate")
  bw  <- get_val("Warming", "estimate")
  bp  <- get_val("Pollution", "estimate")
  bwp <- get_val("Warming:Pollution", "estimate")
  
  pwp <- get_val("Warming:Pollution", "p.value")
  lo  <- get_val("Warming:Pollution", "conf.low")
  hi  <- get_val("Warming:Pollution", "conf.high")
  
  if (any(is.na(c(b0, bw, bp, bwp, pwp, lo, hi)))) {
    stop("extract_barnacle_interaction_class(): Missing fixed-effect term(s) or CI columns in tidy() output.")
  }
  
  d_add <- bw + bp
  d_obs <- d_add + bwp
  
  int_ratio    <- exp(bwp)
  int_ratio_lo <- exp(lo)
  int_ratio_hi <- exp(hi)
  
  overlaps_equiv <- !(int_ratio_hi < equiv_lo || int_ratio_lo > equiv_hi)
  
  dir_label <- dplyr::case_when(
    (d_add > 0  & d_obs < 0) ~ "Reversal (Pos→Neg)",
    (d_add < 0  & d_obs > 0) ~ "Reversal (Neg→Pos)",
    sign(d_add) == sign(bwp) ~ "Synergism (same-sign)",
    TRUE                     ~ "Antagonism (opp-sign)"
  )
  
  class <- dplyr::case_when(
    pwp >= alpha   ~ "Unclassified (p≥0.05)",
    overlaps_equiv ~ "Additive (CI overlaps ±5%)",
    TRUE           ~ dir_label
  )
  
  tibble(
    Date            = label,
    Interaction_Ratio = int_ratio,
    Ratio_CI          = sprintf("[%.3f, %.3f]", int_ratio_lo, int_ratio_hi),
    p_value           = pwp,
    Class             = class,
    Additive_log      = d_add,
    Observed_log      = d_obs,
    Additive_resp     = exp(b0 + d_add),
    Observed_resp     = exp(b0 + d_obs)
  )
}

int_class_tbl <- extract_barnacle_interaction_class(model_final, label = "Whole Season")
print(int_class_tbl)

## 3.4 Conditional plot (two stacked panels; log10 scale with pseudo-count)
# We plot on a log10 axis for readability; pseudo-count is applied only to plotted values.
emm_link_cells <- emmeans(model_final, ~ Warming * Pollution, type = "link")

pred_df <- as.data.frame(emm_link_cells) %>%
  rename(log_pred = emmean, lower_log = lower.CL, upper_log = upper.CL) %>%
  mutate(
    Warming   = factor(Warming,   levels = c(0, 1), labels = c("Ambient", "Warmed")),
    Pollution = factor(Pollution, levels = c(0, 1), labels = c("Non-Polluted", "Polluted")),
    Warming_x = if_else(Warming == "Ambient", 0.4, 0.6)
  )

raw_df <- pad_summary %>%
  mutate(
    Warming   = factor(Warming,   levels = c(0, 1), labels = c("Ambient", "Warmed")),
    Pollution = factor(Pollution, levels = c(0, 1), labels = c("Non-Polluted", "Polluted")),
    Warming_x = if_else(Warming == "Ambient", 0.4, 0.6),
    TreatFill = case_when(
      Pollution == "Non-Polluted" & Warming == "Ambient" ~ "NP_Ambient",
      Pollution == "Non-Polluted" & Warming == "Warmed"  ~ "NP_Warmed",
      Pollution == "Polluted"     & Warming == "Ambient" ~ "P_Ambient",
      Pollution == "Polluted"     & Warming == "Warmed"  ~ "P_Warmed"
    )
  )

site_means <- raw_df %>%
  group_by(Site = `0-Country/site`, Warming, Pollution) %>%
  summarise(site_avg = mean(AvgBarnacles, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Pollution) %>%
  filter(n_distinct(Warming) == 2) %>%
  ungroup() %>%
  mutate(
    Warming_x   = if_else(Warming == "Ambient", 0.4, 0.6),
    site_avg_pc = site_avg + PSEUDO_COUNT,
    TreatFill   = case_when(
      Pollution == "Non-Polluted" & Warming == "Ambient" ~ "NP_Ambient",
      Pollution == "Non-Polluted" & Warming == "Warmed"  ~ "NP_Warmed",
      Pollution == "Polluted"     & Warming == "Ambient" ~ "P_Ambient",
      Pollution == "Polluted"     & Warming == "Warmed"  ~ "P_Warmed"
    )
  )

pred_np <- pred_df %>% filter(Pollution == "Non-Polluted")
pred_p  <- pred_df %>% filter(Pollution == "Polluted")
raw_np  <- raw_df  %>% filter(Pollution == "Non-Polluted")
raw_p   <- raw_df  %>% filter(Pollution == "Polluted")
site_np <- site_means %>% filter(Pollution == "Non-Polluted")
site_p  <- site_means %>% filter(Pollution == "Polluted")

# Back-transform from link:
#   model is log(Avg+1) -> exp(lp) - 1 gives Avg.
# Then add pseudo-count only for plotting on log10 scale.
y_from_link <- function(lp) exp(lp) - 1 + PSEUDO_COUNT

# Non-polluted panel
p_np <- ggplot() +
  geom_ribbon(
    data = pred_np,
    aes(x = Warming_x, ymin = y_from_link(lower_log), ymax = y_from_link(upper_log), group = 1),
    fill = "grey90", alpha = 0.5, colour = NA
  ) +
  # Cross-stratum comparison line (Polluted), faint dashed
  geom_line(
    data = pred_p,
    aes(x = Warming_x, y = y_from_link(log_pred), group = 1),
    colour = "#006400", linewidth = 1.2, linetype = "dashed", alpha = 0.35, lineend = "round"
  ) +
  geom_line(
    data = pred_np,
    aes(x = Warming_x, y = y_from_link(log_pred), group = 1),
    colour = "grey40", linewidth = 2.4, lineend = "round"
  ) +
  geom_jitter(
    data = raw_np,
    aes(x = Warming_x, y = AvgBarnacles + PSEUDO_COUNT, fill = TreatFill),
    width = 0.01, alpha = 0.3, size = 1.8, shape = 21,
    colour = "black", stroke = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = site_np,
    aes(x = Warming_x, y = site_avg_pc, group = Site),
    colour = "grey70", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = site_np,
    aes(x = Warming_x, y = site_avg_pc, fill = TreatFill),
    shape = 21, size = 3.5, stroke = 1.0, colour = "black"
  ) +
  scale_fill_manual(values = TREAT_FILL) +
  scale_x_continuous(
    breaks = c(0.4, 0.6), labels = c("Ambient", "Warmed"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name   = expression(Barnacle~count~"(log"[10]*" scale)"),
    trans  = "log10",
    breaks = c(1, 10, 30, 100, 300),
    labels = label_number(accuracy = 1)
  ) +
  labs(x = "Warming treatment") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 13),
    axis.text       = element_text(size = 11),
    plot.margin     = margin(t = 5, r = 10, b = 2, l = 10)
  )

# Polluted panel
p_p <- ggplot() +
  geom_ribbon(
    data = pred_p,
    aes(x = Warming_x, ymin = y_from_link(lower_log), ymax = y_from_link(upper_log), group = 1),
    fill = "#006400", alpha = 0.25, colour = NA
  ) +
  # Cross-stratum comparison line (Non-Polluted), faint dashed
  geom_line(
    data = pred_np,
    aes(x = Warming_x, y = y_from_link(log_pred), group = 1),
    colour = "grey40", linewidth = 1.2, linetype = "dashed", alpha = 0.35, lineend = "round"
  ) +
  geom_line(
    data = pred_p,
    aes(x = Warming_x, y = y_from_link(log_pred), group = 1),
    colour = "#006400", linewidth = 2.4, lineend = "round"
  ) +
  geom_jitter(
    data = raw_p,
    aes(x = Warming_x, y = AvgBarnacles + PSEUDO_COUNT, fill = TreatFill),
    width = 0.01, alpha = 0.3, size = 1.8, shape = 21,
    colour = "black", stroke = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = site_p,
    aes(x = Warming_x, y = site_avg_pc, group = Site),
    colour = "#006400", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = site_p,
    aes(x = Warming_x, y = site_avg_pc, fill = TreatFill),
    shape = 21, size = 3.5, stroke = 1.0, colour = "black"
  ) +
  scale_fill_manual(values = TREAT_FILL) +
  scale_x_continuous(
    breaks = c(0.4, 0.6), labels = c("Ambient", "Warmed"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name   = expression(Barnacle~count~"(log"[10]*" scale)"),
    trans  = "log10",
    breaks = c(1, 10, 30, 100, 300),
    labels = label_number(accuracy = 1)
  ) +
  labs(x = "Warming treatment") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.y    = element_text(size = 13),
    axis.text       = element_text(size = 11),
    plot.margin     = margin(t = 2, r = 10, b = 5, l = 10)
  )

print(p_np / p_p + plot_layout(heights = c(1, 1)))

## 4. Site-level insights (treatment effects + biogeography) -------------------

# Per-site metadata for ordering (north to south)
site_info <- pad_summary %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site_raw     = `0-Country/site`,
    Site         = clean_site(`0-Country/site`),
    Latitude     = as.numeric(`0-Latitude`),
    Is_temperate = `4-Is_temperate`
  )

# Barnacles support full random slopes: Warming, Pollution, and Interaction
# Extracted directly as BLUP total slopes on the log scale, with SEs derived from the per-site posterior variance-cov matrix.
re_site  <- ranef(model_final, condVar = TRUE)$`0-Country/site`
pv_array <- attr(re_site, "postVar")  # (nTerms x nTerms x nSites)
re_terms <- colnames(re_site)

need_terms <- c("(Intercept)", "Warming", "Pollution", "Warming:Pollution")
if (!all(need_terms %in% re_terms)) {
  stop("Expected random-effect term(s) not found in ranef(model_final)$`0-Country/site`: ",
       paste(setdiff(need_terms, re_terms), collapse = ", "))
}

fix_eff <- fixef(model_final)

# Helper to pull var-cov entry for a given term pair (i, j) across sites
pv_ij <- function(term_i, term_j) {
  ii <- which(re_terms == term_i)
  jj <- which(re_terms == term_j)
  pv_array[ii, jj, ]
}

site_effects <- tibble(
  Site_raw = rownames(re_site),
  Site     = clean_site(rownames(re_site)),
  u_W      = re_site[, "Warming"],
  u_P      = re_site[, "Pollution"],
  u_INT    = re_site[, "Warming:Pollution"]
) %>%
  mutate(
    # Total slopes (log scale)
    warm_log = unname(fix_eff["Warming"]) + u_W,
    poll_log = unname(fix_eff["Pollution"]) + u_P,
    int_log  = unname(fix_eff["Warming:Pollution"]) + u_INT,
    
    # SEs (log scale) from posterior var-cov (per-site vectors)
    warm_se = safe_sqrt(pv_ij("Warming", "Warming")),
    poll_se = safe_sqrt(pv_ij("Pollution", "Pollution")),
    int_se  = safe_sqrt(pv_ij("Warming:Pollution", "Warming:Pollution"))
  )

# Long-form for plotting: one row per site×term
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
    Term      = factor(Term, levels = c("Warming", "Pollution", "Interaction"))
  )

# Global fixed effects
vc_fix <- as.matrix(vcov(model_final))

global_df <- tibble(
  Term       = factor(c("Warming", "Pollution", "Interaction"),
                      levels = c("Warming", "Pollution", "Interaction")),
  slope_log  = c(unname(fix_eff["Warming"]),
                 unname(fix_eff["Pollution"]),
                 unname(fix_eff["Warming:Pollution"])),
  se_log     = c(
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
    Site      = "<b>Global fixed effect</b>"
  )

# Y-order: global at top, then sites north to south
ord_sites <- site_info %>%
  semi_join(all_df %>% distinct(Site), by = "Site") %>%
  arrange(desc(Latitude)) %>%
  pull(Site)

y_levels <- rev(c("<b>Global fixed effect</b>", ord_sites))

# Caterpillar plot: log scale
p_log_barn <- ggplot() +
  geom_rect(
    data  = global_df,
    aes(xmin = lower_log, xmax = upper_log),
    ymin  = -Inf, ymax = Inf,
    fill  = "grey80", alpha = 0.3
  ) +
  geom_vline(
    data = global_df,
    aes(xintercept = slope_log),
    linetype = "dotted",
    colour   = "black"
  ) +
  geom_point(
    data   = global_df,
    aes(x = slope_log, y = Site, fill = Term),
    shape  = 23,
    size   = 5,
    colour = "black"
  ) +
  geom_errorbarh(
    data   = all_df %>% mutate(Site = as.character(Site)),
    aes(y = Site, xmin = lower_log, xmax = upper_log),
    height = 0.2,
    colour = "#666666"
  ) +
  geom_point(
    data = all_df %>% mutate(Site = as.character(Site)),
    aes(x = slope_log, y = Site, colour = Term),
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
    nrow   = 1,
    scales = "free_x",
    labeller = labeller(
      Warming     = "Warming effect",
      Pollution   = "Pollution effect",
      Interaction = "Interaction effect"
    )
  ) +
  scale_x_continuous(name = "Log change in mean barnacle count (log(Avg+1) scale)") +
  scale_y_discrete(
    limits = y_levels,
    labels = pretty_site_labels,
    name   = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text      = element_text(face = "bold"),
    panel.spacing.x = unit(1.5, "lines"),
    axis.text.y     = ggtext::element_markdown(size = 9),
    legend.position = "none"
  )

print(p_log_barn)

## 4.1 Biogeographic tables (latitude + climate zone)

get_site_info <- function(pad_summary) {
  pad_summary %>%
    distinct(`0-Country/site`, .keep_all = TRUE) %>%
    transmute(
      Site_raw = `0-Country/site`,
      Site     = clean_site(`0-Country/site`),
      abs_lat  = abs(as.numeric(`0-Latitude`)),
      
      # Robust to TRUE/FALSE or 1/0 (and keeps NA as NA)
      Climate  = factor(
        dplyr::case_when(
          `4-Is_temperate` %in% c(TRUE, 1)  ~ "Temperate",
          `4-Is_temperate` %in% c(FALSE, 0) ~ "Non-temperate",
          TRUE ~ NA_character_
        ),
        levels = c("Temperate", "Non-temperate")
      )
    )
}

# Pearson correlation with 95% CI; returns NA if too few points/variation
tidy_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; N <- length(x)
  if (N < 3 || dplyr::n_distinct(x) < 2 || dplyr::n_distinct(y) < 2) {
    return(tibble(r = NA_real_, r_L = NA_real_, r_U = NA_real_, p = NA_real_, N = N))
  }
  ct <- cor.test(x, y, method = "pearson", conf.level = 0.95)
  tibble(r = unname(ct$estimate), r_L = ct$conf.int[1], r_U = ct$conf.int[2],
         p = ct$p.value, N = N)
}

# Simple linear slope (y ~ x) with 95% CI; returns NA if underpowered
tidy_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || dplyr::n_distinct(x) < 2) {
    return(tibble(slope = NA_real_, s_L = NA_real_, s_U = NA_real_))
  }
  broom::tidy(lm(y ~ x), conf.int = TRUE) %>%
    dplyr::filter(term == "x") %>%
    transmute(slope = estimate, s_L = conf.low, s_U = conf.high)
}

# Two-sample t-test summary for Temperate vs Non-temperate
tidy_zone <- function(y, zone) {
  d <- tibble(y = y, zone = zone) %>% filter(is.finite(y), !is.na(zone))
  if (dplyr::n_distinct(d$zone) < 2) {
    return(tibble(
      N_Temperate = NA, N_NonTemp = NA, Mean_Temperate = NA, Mean_NonTemp = NA,
      Diff_TminusN = NA, CI_L = NA, CI_U = NA, t = NA, df = NA, p = NA
    ))
  }
  tt <- t.test(y ~ zone, data = d)
  tibble(
    N_Temperate    = sum(d$zone == "Temperate"),
    N_NonTemp      = sum(d$zone == "Non-temperate"),
    Mean_Temperate = mean(d$y[d$zone == "Temperate"]),
    Mean_NonTemp   = mean(d$y[d$zone == "Non-temperate"]),
    Diff_TminusN   = unname(diff(tt$estimate)),
    CI_L = tt$conf.int[1], CI_U = tt$conf.int[2],
    t  = unname(tt$statistic),
    df = unname(tt$parameter),
    p  = unname(tt$p.value)
  )
}

make_biogeog_tables <- function(all_df, pad_summary,
                                response_label = "Barnacles",
                                pct_col = "pct") {
  site_info2 <- get_site_info(pad_summary)
  
  ef <- all_df %>%
    mutate(Site = as.character(Site)) %>%
    filter(Site != "<b>Global fixed effect</b>") %>%
    select(Site, Term, !!pct_col) %>%
    rename(pct_effect = !!pct_col) %>%
    left_join(site_info2, by = "Site") %>%
    filter(is.finite(pct_effect), is.finite(abs_lat), !is.na(Climate))
  
  zone <- ef %>%
    group_by(Term) %>%
    summarise(tidy_zone(pct_effect, Climate), .groups = "drop") %>%
    mutate(Response = response_label) %>%
    relocate(Response, Term) %>%
    mutate(across(c(Mean_Temperate, Mean_NonTemp, Diff_TminusN, t), ~round(., 2)),
           p = signif(p, 3))
  
  latitude <- ef %>%
    group_by(Term) %>%
    summarise(
      tidy_cor(abs_lat, pct_effect),
      tidy_slope(abs_lat, pct_effect),
      .groups = "drop"
    ) %>%
    mutate(Response = response_label) %>%
    relocate(Response, Term) %>%
    mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
           p = signif(p, 3))
  
  list(zone = zone, latitude = latitude)
}

barn_biog <- make_biogeog_tables(all_df, pad_summary,
                                 response_label = "Barnacles",
                                 pct_col = "pct")
print(barn_biog$zone)
print(barn_biog$latitude)

## 5. Meta-regression: site-level warming sensitivity --------------------------

site_warm_base <- tibble(
  Site_raw  = rownames(re_site),
  Site      = clean_site(rownames(re_site)),
  uW        = re_site[, "Warming"],
  uINT      = re_site[, "Warming:Pollution"],
  
  # Per-site posterior var/cov (vectors length = nSites)
  var_uW    = as.numeric(pv_ij("Warming", "Warming")),
  var_uINT  = as.numeric(pv_ij("Warming:Pollution", "Warming:Pollution")),
  cov_uWINT = as.numeric(pv_ij("Warming", "Warming:Pollution"))
) %>%
  mutate(
    se_np = safe_sqrt(var_uW),
    se_p  = safe_sqrt(var_uW + var_uINT + 2 * cov_uWINT)
  )

site_warm <- site_warm_base %>%
  tidyr::expand_grid(Pollution = c("Non-Polluted", "Polluted")) %>%
  mutate(
    slope_log = (unname(fix_eff["Warming"]) + uW) +
      if_else(
        Pollution == "Polluted",
        unname(fix_eff["Warming:Pollution"]) + uINT,
        0
      ),
    
    se_slope = if_else(Pollution == "Non-Polluted", se_np, se_p),
    warm_pct = 100 * (exp(slope_log) - 1)
  ) %>%
  select(Site, Pollution, slope_log, se_slope, warm_pct)

# ΔWarming (°C) from logger summary (hard-coded; see 'Logger processing' supplementary code for summaries)
pad_max_temps <- tibble(
  Site_raw = rep(
    c("Chile","Ecuador","NZ","Denmark","Greenland",
      "Norway","Portugal1","Portugal2","Vietnam1","UK"),
    each = 2
  ),
  group   = rep(c("Ambient", "Warmed"), times = 10),
  MeanMax = c(
    35.40000, 36.82727,
    41.26667, 41.73333,
    36.58182, 36.72500,
    43.13333, 45.75833,
    21.37500, 22.71818,
    31.31818, 34.05833,
    49.56000, 51.04000,
    35.32500, 36.35833,
    50.37500, 51.31000,
    36.44167, 38.05572
  )
) %>%
  mutate(Site = clean_site(Site_raw)) %>%
  tidyr::pivot_wider(names_from = group, values_from = MeanMax) %>%
  transmute(
    Site,
    DeltaWarming = Warmed - Ambient
  )

# Thermal/geography moderators
df_mods <- data_prep %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site    = clean_site(`0-Country/site`),
    abs_lat = abs(`0-Latitude`),
    
    Climate = factor(
      dplyr::case_when(
        `4-Is_temperate` %in% c(TRUE, 1)  ~ "Temperate",
        `4-Is_temperate` %in% c(FALSE, 0) ~ "Non-temperate",
        TRUE ~ NA_character_
      ),
      levels = c("Temperate", "Non-temperate")
    ),
    
    mean_daily_max_temp = `4-mean_daily_max_temp`,
    temp_range_sd       = `4-sd_temp_range`
  )

df_meta_warm <- site_warm %>%
  left_join(pad_max_temps, by = "Site") %>%
  left_join(df_mods, by = "Site")

# 5.1 Warming meta-regression outputs

# Clean subsets
df_dw <- df_meta_warm %>% filter(is.finite(DeltaWarming), is.finite(warm_pct), !is.na(Pollution))
df_np <- df_meta_warm %>% filter(Pollution == "Non-Polluted", is.finite(warm_pct))

df_zone <- df_np %>% filter(!is.na(Climate))
df_lat  <- df_np %>% filter(is.finite(abs_lat))
df_max  <- df_np %>% filter(is.finite(mean_daily_max_temp))
df_sd   <- df_np %>% filter(is.finite(temp_range_sd))

# Table A: warm_pct vs DeltaWarming by stratum
Barnacles_Warm_TableA <- dplyr::bind_rows(
  {
    s <- df_dw %>% filter(Pollution == "Non-Polluted")
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Stratum = "Non-polluted"),
      tidy_cor(s$DeltaWarming, s$warm_pct),
      tidy_slope(s$DeltaWarming, s$warm_pct)
    )
  },
  {
    s <- df_dw %>% filter(Pollution == "Polluted")
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Stratum = "Polluted"),
      tidy_cor(s$DeltaWarming, s$warm_pct),
      tidy_slope(s$DeltaWarming, s$warm_pct)
    )
  }
) %>%
  mutate(
    across(c(r, r_L, r_U, slope, s_L, s_U), ~ round(., 3)),
    p = signif(p, 3)
  )

# Table B: NP-only associations with moderators (lat + thermal metrics)
Barnacles_Warm_TableB <- dplyr::bind_rows(
  {
    s <- df_lat
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Moderator = "Absolute latitude (°)"),
      tidy_cor(s$abs_lat, s$warm_pct),
      tidy_slope(s$abs_lat, s$warm_pct)
    )
  },
  {
    s <- df_max
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Moderator = "Mean daily max (°C)"),
      tidy_cor(s$mean_daily_max_temp, s$warm_pct),
      tidy_slope(s$mean_daily_max_temp, s$warm_pct)
    )
  },
  {
    s <- df_sd
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Moderator = "SD of daily temp range (°C)"),
      tidy_cor(s$temp_range_sd, s$warm_pct),
      tidy_slope(s$temp_range_sd, s$warm_pct)
    )
  }
) %>%
  mutate(
    across(c(r, r_L, r_U, slope, s_L, s_U), ~ round(., 3)),
    p = signif(p, 3)
  )

# Table C: NP-only climate-zone test
Barnacles_Warm_TableC <- dplyr::bind_cols(
  tibble::tibble(Response = "Barnacles", Stratum = "Non-polluted"),
  tidy_zone(df_zone$warm_pct, df_zone$Climate)
) %>%
  relocate(Response, Stratum) %>%
  mutate(
    across(c(Mean_Temperate, Mean_NonTemp, Diff_TminusN, CI_L, CI_U, t), ~ round(., 2)),
    p = signif(p, 3)
  )

print(Barnacles_Warm_TableA)
print(Barnacles_Warm_TableB)
print(Barnacles_Warm_TableC)

# Warming figures
base_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "plain"),
        legend.position = "none")

g_w_dw <- ggplot(df_dw, aes(DeltaWarming, warm_pct)) +
  geom_smooth(aes(color = Pollution, fill = Pollution),
              method = "lm", formula = y ~ x, se = TRUE, alpha = 0.25) +
  geom_point(aes(fill = Pollution),
             shape = 21, size = 2.8, stroke = 1.4, colour = "black") +
  scale_color_manual(values = PAL_STRATA) +
  scale_fill_manual(values  = PAL_STRATA) +
  labs(title = "Warming effect vs ΔWarming (by stratum)",
       x = expression(Delta * "Warming (°C)"),
       y = "% change under warming") +
  base_theme

g_w_zone <- ggplot(df_zone, aes(Climate, warm_pct)) +
  geom_boxplot(width = 0.45, fill = NA, colour = "black", outlier.shape = NA) +
  geom_jitter(
    width = 0.05, shape = 21, size = 2.8, stroke = 1.4,
    fill = "grey75", colour = "black"
  ) +
  labs(title = "Climate zone comparison", x = NULL, y = "% change under warming") +
  base_theme

g_w_lat <- ggplot(df_lat, aes(abs_lat, warm_pct)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4, fill = "grey75", colour = "black") +
  labs(title = "Absolute latitude", x = "Absolute latitude (°)", y = "% change under warming") +
  base_theme

g_w_max <- ggplot(df_max, aes(mean_daily_max_temp, warm_pct)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4, fill = "grey75", colour = "black") +
  labs(title = "Mean daily max", x = "Mean daily max (°C)", y = "% change under warming") +
  base_theme

g_w_sd <- ggplot(df_sd, aes(temp_range_sd, warm_pct)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4, fill = "grey75", colour = "black") +
  labs(title = "SD of daily temp range", x = "SD of daily temp range (°C)", y = "% change under warming") +
  base_theme

g_w_grid <- (g_w_dw | g_w_zone) / (g_w_lat | g_w_max) / (g_w_sd | plot_spacer())
print(g_w_grid)

## 6. Meta-regression: site-level pollution sensitivity ------------------------

site_poll_base <- tibble(
  Site_raw  = rownames(re_site),
  Site      = clean_site(rownames(re_site)),
  uP        = re_site[, "Pollution"],
  uINT      = re_site[, "Warming:Pollution"],
  
  var_uP    = as.numeric(pv_ij("Pollution", "Pollution")),
  var_uINT  = as.numeric(pv_ij("Warming:Pollution", "Warming:Pollution")),
  cov_uPINT = as.numeric(pv_ij("Pollution", "Warming:Pollution"))
) %>%
  mutate(
    se_amb = safe_sqrt(var_uP),
    se_w   = safe_sqrt(var_uP + var_uINT + 2 * cov_uPINT)
  )

site_poll <- site_poll_base %>%
  tidyr::expand_grid(Warming = c("Ambient", "Warmed")) %>%
  mutate(
    slope_log = (unname(fix_eff["Pollution"]) + uP) +
      if_else(
        Warming == "Warmed",
        unname(fix_eff["Warming:Pollution"]) + uINT,
        0
      ),
    
    se_slope = if_else(Warming == "Ambient", se_amb, se_w),
    poll_pct = 100 * (exp(slope_log) - 1)
  ) %>%
  select(Site, Warming, slope_log, se_slope, poll_pct)

# Nitrate enrichment dose as lnR_NO3 (hard-coded from nutrient processing; see 'Pollution processing' supplementary code for summaries).
df_no3_raw <- tibble::tribble(
  ~Site,           ~prop_increase,
  "Chile",          1.0700000,
  "Denmark",        1.2343662,
  "Ecuador",        1.7889979,
  "Greenland",      0.8395716,
  "New Zealand",    2.1634769,
  "Norway",         3.3617540,
  "Portugal 1",    44.0322001,
  "United Kingdom", 1.0363179,
  "Vietnam",        4.0791139
) %>%
  mutate(
    Site    = clean_site(Site),
    lnR_NO3 = log1p(prop_increase)
  )

df_mod_poll <- data_prep %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site      = clean_site(`0-Country/site`),
    human     = `4-pct_human_landuse`,
    precip_mu = `4-mean_daily_precip`
  )

df_meta_poll <- site_poll %>%
  left_join(df_no3_raw %>% select(Site, lnR_NO3), by = "Site") %>%
  left_join(df_mod_poll, by = "Site") %>%
  # Drop Portugal 2 once here
  filter(Site != "Portugal 2")

# lnR meta-regression frame (both warming strata)
df_meta_poll_lnR <- df_meta_poll %>%
  filter(is.finite(lnR_NO3), is.finite(poll_pct), !is.na(Warming))

# Ambient-only moderators frame (human + precip)
df_amb <- df_meta_poll %>%
  filter(Warming == "Ambient") %>%
  filter(is.finite(human), is.finite(precip_mu), is.finite(poll_pct))

## 6.1 Pollution meta-regression outputs ----------------

# Table A: lnR_NO3 vs poll_pct by Warming stratum
Barnacles_Poll_TableA <- dplyr::bind_rows(
  {
    s <- df_meta_poll_lnR %>% filter(Warming == "Ambient")
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Stratum = "Ambient"),
      tidy_cor(s$lnR_NO3, s$poll_pct),
      tidy_slope(s$lnR_NO3, s$poll_pct)
    )
  },
  {
    s <- df_meta_poll_lnR %>% filter(Warming == "Warmed")
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Stratum = "Warmed"),
      tidy_cor(s$lnR_NO3, s$poll_pct),
      tidy_slope(s$lnR_NO3, s$poll_pct)
    )
  }
) %>%
  mutate(
    across(c(r, r_L, r_U, slope, s_L, s_U), ~ round(., 3)),
    p = signif(p, 3)
  )

# Table B: Ambient-only associations with human + precip
Barnacles_Poll_TableB <- dplyr::bind_rows(
  {
    s <- df_amb
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Moderator = "Human land cover (%)"),
      tidy_cor(s$human, s$poll_pct),
      tidy_slope(s$human, s$poll_pct)
    )
  },
  {
    s <- df_amb
    dplyr::bind_cols(
      tibble::tibble(Response = "Barnacles", Moderator = "Mean precipitation (mm/day)"),
      tidy_cor(s$precip_mu, s$poll_pct),
      tidy_slope(s$precip_mu, s$poll_pct)
    )
  }
) %>%
  mutate(
    across(c(r, r_L, r_U, slope, s_L, s_U), ~ round(., 3)),
    p = signif(p, 3)
  )

print(Barnacles_Poll_TableA)
print(Barnacles_Poll_TableB)

## 6.2 Pollution figures -------------------------------------------------------

# Uses base_theme from Section 5 (kept as-is)
g_p_lnR <- ggplot(df_meta_poll_lnR, aes(x = lnR_NO3, y = poll_pct)) +
  geom_smooth(aes(color = Warming, fill = Warming),
              method = "lm", formula = y ~ x, se = TRUE, fullrange = TRUE, alpha = 0.25) +
  geom_point(aes(fill = Warming),
             shape = 21, size = 2.8, stroke = 1.4, colour = "black") +
  scale_color_manual(values = PAL_WARM) +
  scale_fill_manual(values  = PAL_WARM) +
  labs(title = "Pollution effect vs nitrogen response ratio",
       x = expression(ln * "R"[NO[3]] * " (Polluted/Ambient)"),
       y = "% change under pollution") +
  base_theme

g_p_human <- ggplot(df_amb, aes(human, poll_pct)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25, fullrange = TRUE) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Human land cover",
       x = "Human land cover (5 km buffer, %)",
       y = "% change under pollution") +
  base_theme

g_p_prec <- ggplot(df_amb, aes(precip_mu, poll_pct)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25, fullrange = TRUE) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Mean precipitation",
       x = "Mean daily precipitation (mm/day)",
       y = "% change under pollution") +
  base_theme

g_p_grid <- (g_p_lnR | g_p_human) / (g_p_prec | plot_spacer())
print(g_p_grid)

### 6.3 Diagnostics & influence: mean precipitation (Ambient stratum) ----------

# Leave-one-out loop for x and y; returns one row per omitted site
loo_analyse <- function(df, x, y){
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))){
    d <- df[-i, , drop = FALSE]
    rc <- tidy_cor(d[[x]], d[[y]])
    sl <- tidy_slope(d[[x]], d[[y]])
    out[[i]] <- dplyr::bind_cols(
      tibble::tibble(omit_site = df$Site[i], N = rc$N),
      dplyr::select(rc, r, r_L, r_U, p),
      sl
    )
  }
  dplyr::bind_rows(out)
}

# Influence summary for a simple lm used in meta-regression plots
infl_tbl <- function(model, df, label){
  n <- nrow(df); p <- length(stats::coef(model))
  tibble::tibble(
    Site     = df$Site,
    leverage = stats::hatvalues(model),
    cooks_d  = stats::cooks.distance(model),
    rstudent = stats::rstudent(model)
  ) |>
    dplyr::mutate(
      flag_cook = cooks_d > 4/n,         # Cook’s D threshold: > 4/n
      flag_hat  = leverage > (2*p)/n,    # hatvalues threshold: > 2p/n
      Model     = label,
      .before   = 1
    )
}

resp_label <- "Barnacles"

dat_prec <- df_amb |>
  dplyr::select(Site, poll_pct, precip_mu) |>
  dplyr::filter(is.finite(poll_pct), is.finite(precip_mu))

# Baseline correlation/slope
Barnacles_Base_precip_mu <- dplyr::bind_cols(
  tibble::tibble(Response = resp_label,
                 Moderator = "Mean precipitation (mm/day)"),
  tidy_cor(dat_prec$precip_mu, dat_prec$poll_pct),
  tidy_slope(dat_prec$precip_mu, dat_prec$poll_pct)
)

# Leave-one-out (mean precipitation)
Barnacles_LOO_precip_mu <- loo_analyse(dat_prec, "precip_mu", "poll_pct") |>
  dplyr::mutate(Moderator = "Mean precipitation (mm/day)", .before = 1)

# Drop-Vietnam check (labels vary; harmonise using intersect)
vn_hit <- intersect(c("Vietnam 1", "Vietnam"), dat_prec$Site)

if (length(vn_hit) > 0) {
  dat_no_vn <- dplyr::filter(dat_prec, !(Site %in% vn_hit))
  
  Barnacles_DropVietnam <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(Set = "All sites",
                     Moderator = "Mean precipitation (mm/day)"),
      tidy_cor(dat_prec$precip_mu, dat_prec$poll_pct),
      tidy_slope(dat_prec$precip_mu, dat_prec$poll_pct)
    ),
    dplyr::bind_cols(
      tibble::tibble(Set = paste0("No ", paste(vn_hit, collapse = "/")),
                     Moderator = "Mean precipitation (mm/day)"),
      tidy_cor(dat_no_vn$precip_mu, dat_no_vn$poll_pct),
      tidy_slope(dat_no_vn$precip_mu, dat_no_vn$poll_pct)
    )
  ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3)))
} else {
  Barnacles_DropVietnam <- tibble::tibble(
    note = "No Vietnam site label found in Ambient data; no drop-Vietnam summary produced."
  )
}

# Influence diagnostics for single-moderator model
lin_mu_b <- stats::lm(poll_pct ~ precip_mu, data = dat_prec)
Barnacles_Influence_mu <- infl_tbl(lin_mu_b, dat_prec, "poll_pct ~ precip_mu")

Barnacles_Base_precip_mu <- Barnacles_Base_precip_mu |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3)))

Barnacles_LOO_precip_mu <- Barnacles_LOO_precip_mu |>
  dplyr::mutate(dplyr::across(c(r, r_L, r_U, p, slope, s_L, s_U), ~ round(., 3)))

Barnacles_Influence_mu <- Barnacles_Influence_mu |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3)))

print(Barnacles_Base_precip_mu)
print(Barnacles_LOO_precip_mu)
print(Barnacles_DropVietnam)
print(Barnacles_Influence_mu)