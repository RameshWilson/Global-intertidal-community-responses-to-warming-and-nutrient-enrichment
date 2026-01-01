###############################################################################
# Macroalgal cover (Main analyses)
#
#   Plate-level summer means (AvgAlgae = % cover)
#   Final Beta–logit GLMM on proportion (0, 1)
#   * Random effects: (1 + Warming + Pollution + Warming:Pollution | Site)
#   Fixed effects + EMM contrasts
#   Interaction classification
#   Conditional plot
#   Site-level stressor responses (log-odds caterpillar BLUP plot + biogeography tables)
#   Meta-regressions (warming sensitivity; pollution sensitivity)
#   Diagnostics & influence key pollution meta-regression
###############################################################################

## 0. Setup -------------------------------------------------------------------
library(here)
library(tidyverse)
library(glmmTMB)
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
# For beta–logit, exp(beta_int) is an odds ratio (OR).
EQUIV_LO <- 0.95
EQUIV_HI <- 1.05

# Beta-family requires response strictly within (0, 1).
# We clip the proportion using an epsilon so the model is well-defined.
BETA_EPS <- 0.001

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
    stop("standardize_ci(): CI columns not found. Present: ",
         paste(names(df), collapse = ", "))
  }
  
  df %>% dplyr::rename(lower = !!lower_name, upper = !!upper_name)
}

# Variance-covariance arrays can contain tiny negative values due to floating point rounding.
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

## 1.2 Recode treatments
# Warming: Black pad = warmed (1), White pad = ambient (0)
# Pollution: P region = polluted (1), N region = non-polluted (0)
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
    AvgAlgae         = mean(`3-Macroalgae total pct cover`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(AvgAlgae))

## 1.4 Drop sites with extremely low algae signal across summer
# Beta–logit + site-level slopes become unstable when the response is effectively ~0 at a whole site. 
# Such sites contribute little information about proportional changes in cover.
low_sites <- pad_summary %>%
  group_by(`0-Country/site`) %>%
  summarise(mu = mean(AvgAlgae, na.rm = TRUE), .groups = "drop") %>%
  filter(mu < 1) %>%  # <1% average cover across summer
  pull(`0-Country/site`)
print(low_sites) # Chile

## 1.5 Build model response + plot keys
pad_summary <- pad_summary %>%
  filter(!(`0-Country/site` %in% low_sites)) %>%
  mutate(
    # Convert percent cover to proportion
    prop_algae = AvgAlgae / 100,
    
    # Clip to (0, 1) for beta-family compatibility
    prop_algae = pmin(pmax(prop_algae, BETA_EPS), 1 - BETA_EPS),
    
    # Clean site name for joins/plots only (model uses raw site string)
    Site_clean = clean_site(`0-Country/site`)
  )

## 2. Final model specification -----------------------------------------------

## 2.1 Final beta–logit GLMM on prop_algae, with full random effects structure
# Optimiser notes:
#   We set optimizer=optim(method="BFGS") and profile=TRUE to improve stability for the beta–logit + random-slope structure.
algae_beta <- glmmTMB(
  prop_algae ~ Warming * Pollution +
    (1 + Warming + Pollution + Warming:Pollution | `0-Country/site`),
  data   = pad_summary,
  family = beta_family(link = "logit"),
  control = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS"),
    profile   = TRUE
  )
)
summary(algae_beta)

## 2.2 Quick model checks
# Fuller diagnostics / alternative specifications provided in ADDITIONAL.
pdHess_ok <- tryCatch(isTRUE(algae_beta$sdr$pdHess), error = function(e) NA)
conv_code <- tryCatch(algae_beta$fit$convergence,     error = function(e) NA)

model_check_tbl <- tibble(
  pdHess      = pdHess_ok,
  convergence = conv_code
)
print(model_check_tbl)
# Positive-definite Hessian; well-behaved curvature
# Convergence code = 0; successful convergence

## 3. Main results ------------------------------------------------------------

## 3.1 Fixed effects table (logit scale + OR + % odds change)
fe_tbl <- broom.mixed::tidy(algae_beta, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    term,
    estimate_logit  = estimate,
    std_error       = std.error,
    z_value         = statistic,
    p_value         = p.value,
    odds_ratio      = exp(estimate),
    OR_LCL          = exp(conf.low),
    OR_UCL          = exp(conf.high),
    pct_change_odds = 100 * (odds_ratio - 1)
  )
print(fe_tbl)

## 3.2 EMM contrasts
emm_link <- emmeans(algae_beta, ~ Warming * Pollution, type = "link")

# Warmed vs Ambient within each Pollution stratum
warm_contr <- contrast(emm_link, method = "revpairwise", by = "Pollution") %>%
  summary(infer = TRUE) %>%
  standardize_ci() %>%
  transmute(
    Contrast = "Warmed vs Ambient",
    Stratum  = ifelse(Pollution == 0, "Non-Polluted", "Polluted"),
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
    Stratum  = ifelse(Warming == 0, "Ambient", "Warmed"),
    Ratio    = exp(estimate),
    LCL      = exp(lower),
    UCL      = exp(upper),
    `Δ %`    = 100 * (Ratio - 1),
    p        = p.value
  )

emm_table <- bind_rows(warm_contr, poll_contr)
print(emm_table)

## 3.3 Interaction classification (±5% OR band + p-value gate)

# We classify the Warming×Pollution interaction on the model’s natural (logit) scale, but
# report exp(beta_int) as an odds ratio (OR) using:
#   - p<0.05 gate on the interaction term
#   - ±5% equivalence band on OR: [0.95, 1.05] to determine meaningful deviations from null
# Directional labels use additive vs observed changes on the logit scale:
#   d_add = beta_W + beta_P
#   d_obs = d_add + beta_WP
extract_interaction_class_logit <- function(model,
                                            alpha = 0.05,
                                            equiv_lo = EQUIV_LO,
                                            equiv_hi = EQUIV_HI,
                                            label = "Whole Season") {
  
  fe <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  getv <- function(term, col) {
    out <- fe %>% filter(.data$term == term) %>% slice(1) %>% pull(.data[[col]])
    if (length(out) == 0) NA_real_ else out
  }
  
  b0  <- getv("(Intercept)", "estimate")
  bw  <- getv("Warming", "estimate")
  bp  <- getv("Pollution", "estimate")
  bwp <- getv("Warming:Pollution", "estimate")
  
  pwp <- getv("Warming:Pollution", "p.value")
  lo  <- getv("Warming:Pollution", "conf.low")
  hi  <- getv("Warming:Pollution", "conf.high")
  
  if (any(is.na(c(b0, bw, bp, bwp, pwp, lo, hi)))) {
    stop("extract_interaction_class_logit(): Missing fixed-effect term(s) or CI columns.")
  }
  
  d_add <- bw + bp
  d_obs <- d_add + bwp
  
  OR    <- exp(bwp)
  OR_lo <- exp(lo)
  OR_hi <- exp(hi)
  
  overlaps_equiv <- !(OR_hi < equiv_lo || OR_lo > equiv_hi)
  
  dir_label <- dplyr::case_when(
    (d_add > 0 & d_obs < 0) ~ "Reversal (Pos to Neg)",
    (d_add < 0 & d_obs > 0) ~ "Reversal (Neg to Pos)",
    sign(d_add) == sign(bwp) ~ "Synergism (same-sign)",
    TRUE ~ "Antagonism (opp-sign)"
  )
  
  class <- dplyr::case_when(
    is.na(pwp) | pwp >= alpha ~ "Unclassified (p≥0.05)",
    overlaps_equiv            ~ "Additive (CI overlaps ±5%)",
    TRUE                      ~ dir_label
  )
  
  tibble(
    Date = label,
    Interaction_OR = OR,
    OR_CI = sprintf("[%.3f, %.3f]", OR_lo, OR_hi),
    p_value = pwp,
    Class = class,
    Additive_logit = d_add,
    Observed_logit = d_obs,
    Additive_prob  = plogis(b0 + d_add),
    Observed_prob  = plogis(b0 + d_obs)
  )
}

int_class_tbl <- extract_interaction_class_logit(algae_beta, label = "Whole Season")
print(int_class_tbl)

## 3.4 Conditional plot
emm_resp <- emmeans(algae_beta, ~ Warming * Pollution, type = "response")

pred <- as.data.frame(emm_resp) %>%
  standardize_ci() %>%
  rename(prop_pred = response, lower_prop = lower, upper_prop = upper) %>%
  mutate(
    Warming     = factor(Warming,   levels = c(0, 1), labels = c("Ambient", "Warmed")),
    Pollution   = factor(Pollution, levels = c(0, 1), labels = c("Non-Polluted", "Polluted")),
    Warming_num = if_else(Warming == "Ambient", 0.4, 0.6)
  )

raw_df <- pad_summary %>%
  mutate(
    Warming     = factor(Warming,   levels = c(0, 1), labels = c("Ambient", "Warmed")),
    Pollution   = factor(Pollution, levels = c(0, 1), labels = c("Non-Polluted", "Polluted")),
    Warming_num = if_else(Warming == "Ambient", 0.4, 0.6),
    TreatFill   = case_when(
      Pollution == "Non-Polluted" & Warming == "Ambient" ~ "NP_Ambient",
      Pollution == "Non-Polluted" & Warming == "Warmed"  ~ "NP_Warmed",
      Pollution == "Polluted"     & Warming == "Ambient" ~ "P_Ambient",
      Pollution == "Polluted"     & Warming == "Warmed"  ~ "P_Warmed"
    )
  )

site_means <- raw_df %>%
  group_by(Site = `0-Country/site`, Warming, Pollution) %>%
  summarise(site_prop = mean(prop_algae, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site, Pollution) %>%
  filter(n_distinct(Warming) == 2) %>%
  ungroup() %>%
  mutate(
    Warming_num = if_else(Warming == "Ambient", 0.4, 0.6),
    TreatFill   = case_when(
      Pollution == "Non-Polluted" & Warming == "Ambient" ~ "NP_Ambient",
      Pollution == "Non-Polluted" & Warming == "Warmed"  ~ "NP_Warmed",
      Pollution == "Polluted"     & Warming == "Ambient" ~ "P_Ambient",
      Pollution == "Polluted"     & Warming == "Warmed"  ~ "P_Warmed"
    )
  )

pred_np <- pred %>% filter(Pollution == "Non-Polluted")
pred_p  <- pred %>% filter(Pollution == "Polluted")

raw_np  <- raw_df %>% filter(Pollution == "Non-Polluted")
raw_p   <- raw_df %>% filter(Pollution == "Polluted")

site_np <- site_means %>% filter(Pollution == "Non-Polluted")
site_p  <- site_means %>% filter(Pollution == "Polluted")

# Non-polluted panel
p_np <- ggplot() +
  geom_ribbon(
    data = pred_np,
    aes(x = Warming_num, ymin = lower_prop, ymax = upper_prop),
    fill = "grey90", alpha = 0.5, colour = NA
  ) +
  # Cross-stratum comparison line (Polluted), faint dashed
  geom_line(
    data = pred_p,
    aes(x = Warming_num, y = prop_pred, group = 1),
    colour = "#004d00", linewidth = 1.2, linetype = "dashed", alpha = 0.35, lineend = "round"
  ) +
  geom_line(
    data = pred_np,
    aes(x = Warming_num, y = prop_pred, group = 1),
    colour = "grey40", linewidth = 2.4, lineend = "round"
  ) +
  geom_jitter(
    data = raw_np,
    aes(x = Warming_num, y = prop_algae, fill = TreatFill),
    width = 0.03, alpha = 0.4, size = 1.8, shape = 21,
    colour = "black", stroke = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = site_np,
    aes(x = Warming_num, y = site_prop, group = Site),
    colour = "grey70", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = site_np,
    aes(x = Warming_num, y = site_prop, fill = TreatFill),
    shape = 21, size = 3.5, stroke = 1.0, colour = "black"
  ) +
  scale_fill_manual(values = TREAT_FILL) +
  scale_x_continuous(
    breaks = c(0.4, 0.6), labels = c("Ambient", "Warmed"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name   = "Macroalgal cover (%)",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent
  ) +
  labs(x = "Warming treatment") +
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
    aes(x = Warming_num, ymin = lower_prop, ymax = upper_prop),
    fill = "#006400", alpha = 0.25, colour = NA
  ) +
  # Cross-stratum comparison line (Non-Polluted), faint dashed
  geom_line(
    data = pred_np,
    aes(x = Warming_num, y = prop_pred, group = 1),
    colour = "grey40", linewidth = 1.2, linetype = "dashed", alpha = 0.35, lineend = "round"
  ) +
  geom_line(
    data = pred_p,
    aes(x = Warming_num, y = prop_pred, group = 1),
    colour = "#004d00", linewidth = 2.4, lineend = "round"
  ) +
  geom_jitter(
    data = raw_p,
    aes(x = Warming_num, y = prop_algae, fill = TreatFill),
    width = 0.03, alpha = 0.4, size = 1.8, shape = 21,
    colour = "black", stroke = 0.3, show.legend = FALSE
  ) +
  geom_line(
    data = site_p,
    aes(x = Warming_num, y = site_prop, group = Site),
    colour = "#006400", linewidth = 0.4, alpha = 0.6
  ) +
  geom_point(
    data = site_p,
    aes(x = Warming_num, y = site_prop, fill = TreatFill),
    shape = 21, size = 3.5, stroke = 1.0, colour = "black"
  ) +
  scale_fill_manual(values = TREAT_FILL) +
  scale_x_continuous(
    breaks = c(0.4, 0.6), labels = c("Ambient", "Warmed"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name   = "Macroalgal cover (%)",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent
  ) +
  labs(x = "Warming treatment") +
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

# Random effects (BLUPs) and their conditional variance-covariance arrays
re_site  <- ranef(algae_beta, condVar = TRUE)$cond$`0-Country/site`
cv_array <- attr(re_site, "condVar")     # [nTerms x nTerms x nSites]
re_terms <- colnames(re_site)

need_terms <- c("(Intercept)", "Warming", "Pollution", "Warming:Pollution")
if (!all(need_terms %in% re_terms)) {
  stop("Expected random-effect term(s) not found in ranef(algae_beta)$cond$`0-Country/site`: ",
       paste(setdiff(need_terms, re_terms), collapse = ", "))
}

# Fixed effects and their vcov matrix (conditional component)
fix_eff <- fixef(algae_beta)$cond
vc_fix  <- vcov(algae_beta)$cond

# Helper to pull a variance/covariance entry (term_i, term_j) across sites
pv_ij <- function(term_i, term_j) {
  ii <- which(re_terms == term_i)
  jj <- which(re_terms == term_j)
  cv_array[ii, jj, ]
}

site_effects <- tibble(
  Site_raw = rownames(re_site),
  Site     = clean_site(rownames(re_site)),
  u_W      = re_site[, "Warming"],
  u_P      = re_site[, "Pollution"],
  u_INT    = re_site[, "Warming:Pollution"]
) %>%
  mutate(
    # Total slopes (logit / log-odds scale)
    warm_logit = unname(fix_eff["Warming"]) + u_W,
    poll_logit = unname(fix_eff["Pollution"]) + u_P,
    int_logit  = unname(fix_eff["Warming:Pollution"]) + u_INT,
    
    # SEs for each slope (from conditional var-cov diagonals)
    warm_se = safe_sqrt(pv_ij("Warming", "Warming")),
    poll_se = safe_sqrt(pv_ij("Pollution", "Pollution")),
    int_se  = safe_sqrt(pv_ij("Warming:Pollution", "Warming:Pollution"))
  )

# Long-form for plotting: one row per site×term
all_df <- bind_rows(
  site_effects %>% transmute(Site, Term = "Warming",     slope_link = warm_logit, se_link = warm_se),
  site_effects %>% transmute(Site, Term = "Pollution",   slope_link = poll_logit, se_link = poll_se),
  site_effects %>% transmute(Site, Term = "Interaction", slope_link = int_logit,  se_link = int_se)
) %>%
  mutate(
    lower_link = slope_link - 1.96 * se_link,
    upper_link = slope_link + 1.96 * se_link,
    
    # For interpretability: % change in odds corresponding to a log-odds slope
    pct_odds   = 100 * (exp(slope_link) - 1),
    lower_pct  = 100 * (exp(lower_link) - 1),
    upper_pct  = 100 * (exp(upper_link) - 1),
    Term       = factor(Term, levels = c("Warming", "Pollution", "Interaction"))
  )

# Global fixed effects (for the grey reference band at top of each facet)
global_df <- tibble(
  Term       = factor(c("Warming", "Pollution", "Interaction"),
                      levels = c("Warming", "Pollution", "Interaction")),
  slope_link = c(unname(fix_eff["Warming"]),
                 unname(fix_eff["Pollution"]),
                 unname(fix_eff["Warming:Pollution"])),
  se_link    = c(
    safe_sqrt(as.matrix(vc_fix)["Warming", "Warming"]),
    safe_sqrt(as.matrix(vc_fix)["Pollution", "Pollution"]),
    safe_sqrt(as.matrix(vc_fix)["Warming:Pollution", "Warming:Pollution"])
  )
) %>%
  mutate(
    lower_link = slope_link - 1.96 * se_link,
    upper_link = slope_link + 1.96 * se_link,
    pct_odds   = 100 * (exp(slope_link) - 1),
    lower_pct  = 100 * (exp(lower_link) - 1),
    upper_pct  = 100 * (exp(upper_link) - 1),
    Site       = "<b>Global fixed effect</b>"
  )

# Y-order: global at top, then sites north to south
ord_sites <- site_info %>%
  semi_join(all_df %>% distinct(Site), by = "Site") %>%
  arrange(desc(Latitude)) %>%
  pull(Site)

y_levels <- rev(c("<b>Global fixed effect</b>", ord_sites))

# Caterpillar plot: log-odds scale
p_logit_algae <- ggplot() +
  geom_rect(
    data  = global_df,
    aes(xmin = lower_link, xmax = upper_link),
    ymin  = -Inf, ymax = Inf,
    fill  = "grey80", alpha = 0.3
  ) +
  geom_vline(
    data = global_df,
    aes(xintercept = slope_link),
    linetype = "dotted",
    colour   = "black"
  ) +
  geom_point(
    data   = global_df,
    aes(x = slope_link, y = Site, fill = Term),
    shape  = 23,
    size   = 5,
    colour = "black"
  ) +
  geom_errorbarh(
    data   = all_df %>% mutate(Site = as.character(Site)),
    aes(y = Site, xmin = lower_link, xmax = upper_link),
    height = 0.2,
    colour = "#666666"
  ) +
  geom_point(
    data = all_df %>% mutate(Site = as.character(Site)),
    aes(x = slope_link, y = Site, colour = Term),
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
  scale_x_continuous(name = "Log-odds change in mean macroalgal cover") +
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

print(p_logit_algae)

## 4.1 Biogeographic tables (latitude + climate zone)
get_site_info <- function(pad_summary) {
  pad_summary %>%
    distinct(`0-Country/site`, .keep_all = TRUE) %>%
    transmute(
      Site_raw = `0-Country/site`,
      Site     = clean_site(`0-Country/site`),
      abs_lat  = abs(as.numeric(`0-Latitude`)),
      Climate  = factor(ifelse(`4-Is_temperate`, "Temperate", "Non-temperate"),
                        levels = c("Temperate", "Non-temperate"))
    )
}

# Pearson correlation with 95% CI; returns NA if too few points/variation
tidy_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  N <- length(x)
  if (N < 3 || dplyr::n_distinct(x) < 2 || dplyr::n_distinct(y) < 2) {
    return(tibble(r = NA_real_, r_L = NA_real_, r_U = NA_real_, p = NA_real_, N = N))
  }
  ct <- cor.test(x, y, method = "pearson", conf.level = 0.95)
  tibble(r = unname(ct$estimate), r_L = ct$conf.int[1], r_U = ct$conf.int[2],
         p = unname(ct$p.value), N = N)
}

# Simple linear slope (y ~ x) with 95% CI; returns NA if underpowered
tidy_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || dplyr::n_distinct(x) < 2) {
    return(tibble(slope = NA_real_, s_L = NA_real_, s_U = NA_real_))
  }
  broom::tidy(lm(y ~ x), conf.int = TRUE) %>%
    filter(term == "x") %>%
    transmute(slope = estimate, s_L = conf.low, s_U = conf.high)
}

# Two-sample t-test summary for Temperate vs Non-temperate
tidy_zone <- function(y, zone) {
  d <- tibble(y = y, zone = zone) %>% filter(is.finite(y), !is.na(zone))
  if (dplyr::n_distinct(d$zone) < 2) {
    return(tibble(
      N_Temperate = NA, N_NonTemp = NA,
      Mean_Temperate = NA, Mean_NonTemp = NA,
      Diff_TminusN = NA, CI_L = NA, CI_U = NA,
      t = NA, df = NA, p = NA
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

make_biogeog_tables <- function(effects_df, pad_summary, response_label, pct_col) {
  site_info2 <- get_site_info(pad_summary)
  
  ef <- effects_df %>%
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
    summarise(tidy_cor(abs_lat, pct_effect), tidy_slope(abs_lat, pct_effect), .groups = "drop") %>%
    mutate(Response = response_label) %>%
    relocate(Response, Term) %>%
    mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
           p = signif(p, 3))
  
  list(zone = zone, latitude = latitude)
}

alg_biog <- make_biogeog_tables(
  effects_df = all_df,
  pad_summary = pad_summary,
  response_label = "Macroalgae",
  pct_col = "pct_odds"
)

print(alg_biog$zone)
print(alg_biog$latitude)

## 5. Meta-regression: site-level warming sensitivity --------------------------

re_df <- tibble(
  Site_raw = rownames(re_site),
  Site     = clean_site(rownames(re_site)),
  uW       = re_site[, "Warming"],
  uP       = re_site[, "Pollution"],
  uINT     = re_site[, "Warming:Pollution"]
)

df_warm <- re_df %>%
  tidyr::expand_grid(Pollution = c("Non-Polluted", "Polluted")) %>%
  mutate(
    slope_logit = (unname(fix_eff["Warming"]) + uW) +
      if_else(Pollution == "Polluted",
              unname(fix_eff["Warming:Pollution"]) + uINT, 0),
    warm_pct_odds = 100 * (exp(slope_logit) - 1)
  )

sites_in_model <- tibble(Site = clean_site(rownames(re_site)))

# ΔWarming (°C) from logger summary (hard-coded; see logger-processing code for provenance).
pad_max_temps <- tibble::tribble(
  ~Site,         ~Ambient, ~Warmed,
  "Argentina",    29.55833, 31.08333,
  "Ecuador",      41.26667, 41.73333,
  "New Zealand",  36.58182, 36.72500,
  "Denmark",      43.13333, 45.75833,
  "Greenland",    21.37500, 22.71818,
  "Norway",       31.31818, 34.05833,
  "Portugal 1",   49.56000, 51.04000,
  "Portugal 2",   35.32500, 36.35833,
  "Vietnam",      50.37500, 51.31000,
  "UK",           36.44167, 38.05572
) %>%
  mutate(Site = clean_site(Site), DeltaWarming = Warmed - Ambient) %>%
  select(Site, DeltaWarming) %>%
  semi_join(sites_in_model, by = "Site")

# Thermal/geography moderators
df_mods <- data_prep %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site = clean_site(`0-Country/site`),
    abs_lat = abs(`0-Latitude`),
    Climate = factor(ifelse(`4-Is_temperate`, "Temperate", "Non-temperate"),
                     levels = c("Temperate", "Non-temperate")),
    mean_daily_max_temp = `4-mean_daily_max_temp`,
    temp_range_sd       = `4-sd_temp_range`
  ) %>%
  semi_join(sites_in_model, by = "Site")

df_meta_warm <- df_warm %>%
  left_join(pad_max_temps, by = "Site") %>%
  left_join(df_mods,       by = "Site")

# 5.1 Warming meta-regression outputs

# Table A: ΔWarming association in each Pollution stratum
Algae_Warm_TableA <- bind_rows(
  { s <- df_meta_warm %>% filter(Pollution == "Non-Polluted");
  bind_cols(
    tibble(Response = "Macroalgae", Stratum = "Non-Polluted"),
    tidy_cor(s$DeltaWarming, s$warm_pct_odds),
    tidy_slope(s$DeltaWarming, s$warm_pct_odds)
  )
  },
  { s <- df_meta_warm %>% filter(Pollution == "Polluted");
  bind_cols(
    tibble(Response = "Macroalgae", Stratum = "Polluted"),
    tidy_cor(s$DeltaWarming, s$warm_pct_odds),
    tidy_slope(s$DeltaWarming, s$warm_pct_odds)
  )
  }
) %>%
  mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
         p = signif(p, 3))

# Table B: NP-only associations with moderators (lat + thermal metrics)
df_np <- df_meta_warm %>%
  filter(Pollution == "Non-Polluted", is.finite(warm_pct_odds))

Algae_Warm_TableB <- bind_rows(
  { s <- df_np %>% filter(is.finite(abs_lat));
  bind_cols(
    tibble(Response = "Macroalgae", Moderator = "Absolute latitude (°)"),
    tidy_cor(s$abs_lat, s$warm_pct_odds),
    tidy_slope(s$abs_lat, s$warm_pct_odds)
  )
  },
  { s <- df_np %>% filter(is.finite(mean_daily_max_temp));
  bind_cols(
    tibble(Response = "Macroalgae", Moderator = "Mean daily max (°C)"),
    tidy_cor(s$mean_daily_max_temp, s$warm_pct_odds),
    tidy_slope(s$mean_daily_max_temp, s$warm_pct_odds)
  )
  },
  { s <- df_np %>% filter(is.finite(temp_range_sd));
  bind_cols(
    tibble(Response = "Macroalgae", Moderator = "SD of daily temp range (°C)"),
    tidy_cor(s$temp_range_sd, s$warm_pct_odds),
    tidy_slope(s$temp_range_sd, s$warm_pct_odds)
  )
  }
) %>%
  mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
         p = signif(p, 3))

# Table C: NP-only climate-zone test
df_zone <- df_np %>% filter(!is.na(Climate), is.finite(warm_pct_odds))
Algae_Warm_TableC <- bind_cols(
  tibble(Response = "Macroalgae", Stratum = "Non-Polluted"),
  tidy_zone(df_zone$warm_pct_odds, df_zone$Climate)
) %>%
  relocate(Response, Stratum) %>%
  mutate(across(c(Mean_Temperate, Mean_NonTemp, Diff_TminusN, CI_L, CI_U, t),
                ~round(., 2)),
         p = signif(p, 3))

print(Algae_Warm_TableA)
print(Algae_Warm_TableB)
print(Algae_Warm_TableC)

# Warming figures
df_dw   <- df_meta_warm %>% filter(is.finite(DeltaWarming), is.finite(warm_pct_odds), !is.na(Pollution))
df_lat  <- df_np %>% filter(is.finite(abs_lat))
df_max  <- df_np %>% filter(is.finite(mean_daily_max_temp))
df_sd   <- df_np %>% filter(is.finite(temp_range_sd))
df_zone <- df_np %>% filter(!is.na(Climate))

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "plain"),
        legend.position = "none")

g_w_dw <- ggplot(df_dw, aes(DeltaWarming, warm_pct_odds)) +
  geom_smooth(aes(color = Pollution, fill = Pollution),
              method = "lm", formula = y ~ x, se = TRUE, alpha = 0.25) +
  geom_point(aes(fill = Pollution),
             shape = 21, size = 2.8, stroke = 1.4, colour = "black") +
  scale_color_manual(values = PAL_STRATA) +
  scale_fill_manual(values  = PAL_STRATA) +
  labs(title = "Warming effect vs ΔWarming (odds; by stratum)",
       x = expression(Delta * "Warming (°C)"),
       y = "% odds change under warming") +
  base_theme

g_w_zone <- ggplot(df_zone, aes(Climate, warm_pct_odds)) +
  geom_boxplot(width = 0.45, fill = NA, colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, shape = 21, size = 2.8, stroke = 1.4,
              fill = "grey75", colour = "black") +
  labs(title = "Climate zone comparison (NP only)", x = NULL, y = "% odds change under warming") +
  base_theme

g_w_lat <- ggplot(df_lat, aes(abs_lat, warm_pct_odds)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Absolute latitude (NP only)", x = "Absolute latitude (°)", y = "% odds change under warming") +
  base_theme

g_w_max <- ggplot(df_max, aes(mean_daily_max_temp, warm_pct_odds)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Mean daily max (NP only)", x = "Mean daily max (°C)", y = "% odds change under warming") +
  base_theme

g_w_sd <- ggplot(df_sd, aes(temp_range_sd, warm_pct_odds)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "SD of daily temp range (NP only)", x = "SD of daily temp range (°C)", y = "% odds change under warming") +
  base_theme

print((g_w_dw | g_w_zone) / (g_w_lat | g_w_max) / (g_w_sd | plot_spacer()))

## 6. Meta-regression: site-level pollution sensitivity ------------------------

df_poll <- re_df %>%
  tidyr::expand_grid(Warming = c("Ambient", "Warmed")) %>%
  mutate(
    slope_logit = (unname(fix_eff["Pollution"]) + uP) +
      if_else(Warming == "Warmed",
              unname(fix_eff["Warming:Pollution"]) + uINT, 0),
    poll_pct_odds = 100 * (exp(slope_logit) - 1)
  )

sites_in_model <- tibble(Site = clean_site(rownames(re_site)))

# Nitrate enrichment dose as lnR_NO3 (hard-coded from nutrient processing; see 'Pollution processing' supplementary code for summaries).
df_no3 <- tibble::tribble(
  ~Site,           ~prop_increase,
  "Argentina",      1.8850390,
  "Denmark",        1.2343662,
  "Ecuador",        1.7889979,
  "Greenland",      0.8395716,
  "New Zealand",    2.1634769,
  "Norway",         3.3617540,
  "Portugal 1",    44.0322001,
  "United Kingdom", 1.0363179,
  "Vietnam",        4.0791139
) %>%
  mutate(Site = clean_site(Site), lnR_NO3 = log1p(prop_increase)) %>%
  select(Site, lnR_NO3) %>%
  semi_join(sites_in_model, by = "Site") %>%
  filter(Site != "Portugal 2")

df_mod_poll <- data_prep %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site      = clean_site(`0-Country/site`),
    human     = `4-pct_human_landuse`,
    precip_mu = `4-mean_daily_precip`
  ) %>%
  semi_join(sites_in_model, by = "Site")

df_meta_poll <- df_poll %>%
  left_join(df_no3,      by = "Site") %>%
  left_join(df_mod_poll, by = "Site")

## 6.1 Pollution meta-regression outputs ----------------

# Table A: lnR_NO3 vs poll_pct by Warming stratum
Algae_Poll_TableA <- bind_rows(
  { s <- df_meta_poll %>% filter(Warming == "Ambient", is.finite(lnR_NO3), is.finite(poll_pct_odds));
  bind_cols(
    tibble(Response = "Macroalgae", Stratum = "Ambient"),
    tidy_cor(s$lnR_NO3, s$poll_pct_odds),
    tidy_slope(s$lnR_NO3, s$poll_pct_odds)
  )
  },
  { s <- df_meta_poll %>% filter(Warming == "Warmed", is.finite(lnR_NO3), is.finite(poll_pct_odds));
  bind_cols(
    tibble(Response = "Macroalgae", Stratum = "Warmed"),
    tidy_cor(s$lnR_NO3, s$poll_pct_odds),
    tidy_slope(s$lnR_NO3, s$poll_pct_odds)
  )
  }
) %>%
  mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
         p = signif(p, 3))

# Table B: Ambient-only associations with human + precip
df_amb <- df_meta_poll %>%
  filter(Warming == "Ambient") %>%
  filter(is.finite(human), is.finite(precip_mu), is.finite(poll_pct_odds))

Algae_Poll_TableB <- bind_rows(
  { s <- df_amb %>% filter(is.finite(human));
  bind_cols(
    tibble(Response = "Macroalgae", Moderator = "Human land cover (%)", Stratum = "Ambient"),
    tidy_cor(s$human, s$poll_pct_odds),
    tidy_slope(s$human, s$poll_pct_odds)
  )
  },
  { s <- df_amb %>% filter(is.finite(precip_mu));
  bind_cols(
    tibble(Response = "Macroalgae", Moderator = "Mean precipitation (mm/day)", Stratum = "Ambient"),
    tidy_cor(s$precip_mu, s$poll_pct_odds),
    tidy_slope(s$precip_mu, s$poll_pct_odds)
  )
  }
) %>%
  mutate(across(c(r, r_L, r_U, slope, s_L, s_U), ~round(., 3)),
         p = signif(p, 3))

print(Algae_Poll_TableA)
print(Algae_Poll_TableB)

## 6.2 Pollution figures -------------------------------------------------------

df_meta_poll_lnR <- df_meta_poll %>%
  filter(is.finite(lnR_NO3), is.finite(poll_pct_odds)) %>%
  filter(Site != "Portugal 2")

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "plain"),
        legend.position = "none")

g_p_lnR <- ggplot(df_meta_poll_lnR, aes(lnR_NO3, poll_pct_odds)) +
  geom_smooth(aes(color = Warming, fill = Warming),
              method = "lm", formula = y ~ x, se = TRUE, fullrange = TRUE, alpha = 0.25) +
  geom_point(aes(fill = Warming),
             shape = 21, size = 2.8, stroke = 1.4, colour = "black") +
  scale_color_manual(values = PAL_WARM) +
  scale_fill_manual(values  = PAL_WARM) +
  labs(title = "Pollution effect vs nitrogen response ratio (odds)",
       x = expression(ln * "R"[NO[3]] * " (Polluted/Ambient)"),
       y = "% odds change under pollution") +
  base_theme

g_p_human <- ggplot(df_amb, aes(human, poll_pct_odds)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25, fullrange = TRUE) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Human land cover (Ambient only)",
       x = "Human land cover (5 km buffer, %)",
       y = "% odds change under pollution") +
  base_theme

g_p_prec <- ggplot(df_amb, aes(precip_mu, poll_pct_odds)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", fill = "grey80", alpha = 0.25, fullrange = TRUE) +
  geom_point(shape = 21, size = 2.8, stroke = 1.4,
             fill = "grey75", colour = "black") +
  labs(title = "Mean precipitation (Ambient only)",
       x = "Mean daily precipitation (mm/day)",
       y = "% odds change under pollution") +
  base_theme

print((g_p_lnR | g_p_human) / (g_p_prec | plot_spacer()))