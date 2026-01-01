###############################################################################
# Algae (additional analyses)
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
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(broom.mixed)
library(patchwork)
library(ggnewscale)
library(ggtext)
library(grid)
library(scales)

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
    x == "Vietnam1"         ~ "Vietnam 1",
    x == "Vietnam"          ~ "Vietnam 1",
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

# DHARMa and pearson-style diagnostics and residuals for glmmTBM
make_diag_panel_glmmTMB <- function(mod, data_for_groups, name = "Model") {
  sim <- simulateResiduals(mod, n = 1000)
  
  stats <- c(
    KS   = testUniformity(sim)$p.value,
    DISP = testDispersion(sim)$p.value,
    ZI   = testZeroInflation(sim)$p.value,
    OUT  = testOutliers(sim)$p.value
  )
  
  dfD <- tibble(
    scaled = sim$scaledResiduals,
    grp    = data_for_groups$group
  )
  
  p1 <- ggplot(dfD, aes(sample = scaled)) +
    stat_qq() + stat_qq_line() +
    labs(
      title = paste0(name, " — DHARMa QQ"),
      subtitle = sprintf(
        "KS p=%.3f | Disp p=%.3f | ZI p=%.3f | Out p=%.3f",
        stats["KS"], stats["DISP"], stats["ZI"], stats["OUT"]
      )
    ) +
    theme_minimal()
  
  p2 <- ggplot(dfD, aes(x = grp, y = scaled)) +
    geom_boxplot(fill = "grey90", colour = "black") +
    labs(
      title = paste0(name, " — DHARMa by group"),
      x = "Warming × Pollution",
      y = "Scaled residuals"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dfP <- data_for_groups %>%
    mutate(
      .fitted  = predict(mod, type = "response"),
      .resid   = residuals(mod, type = "pearson"),
      abs_sqrt = sqrt(abs(.resid))
    )
  
  p3 <- ggplot(dfP, aes(sample = .resid)) +
    stat_qq() + stat_qq_line() +
    labs(title = paste0(name, " — Pearson QQ"), y = "Pearson residual") +
    theme_minimal()
  
  p4 <- ggplot(dfP, aes(x = .fitted, y = .resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste0(name, " — Residuals vs fitted"),
      x = "Fitted (response scale)",
      y = "Pearson residual"
    ) +
    theme_minimal()
  
  p5 <- ggplot(dfP, aes(.resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.7) +
    labs(
      title = paste0(name, " — Pearson residual histogram"),
      x = "Pearson residual", y = "Density"
    ) +
    theme_minimal()
  
  p6 <- ggplot(dfP, aes(x = .fitted, y = abs_sqrt)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      title = paste0(name, " — Scale–location"),
      x = "Fitted (response scale)",
      y = "√|Pearson residual|"
    ) +
    theme_minimal()
  
  (p1 | p2) / (p3 | p4) / (p5 | p6) +
    plot_annotation(title = paste("Diagnostics:", name))
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
    AvgAlgae         = mean(`3-Macroalgae total pct cover`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(AvgAlgae))

## 1.4 Drop 'mean < 1%' sites
low_algae_sites <- pad_summary %>%
  group_by(`0-Country/site`) %>%
  summarise(mean_algae_site = mean(AvgAlgae, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_algae_site < 1) %>%
  pull(`0-Country/site`)

## 1.5 Build model response + plot keys
pad_summary <- pad_summary %>%
  filter(!(`0-Country/site` %in% low_algae_sites)) %>%
  mutate(
    # Convert % cover to proportion and clip to (0,1) for beta likelihood.
    prop_algae = AvgAlgae / 100,
    prop_algae = pmin(pmax(prop_algae, 0.001), 0.999),
    
    Site_clean = clean_site(`0-Country/site`),
    group      = interaction(Warming, Pollution, sep = "×"),
    obsID      = row_number()
  )

## 1.6 Exploratory raw visualisation
pad_summary <- add_treatment_factor(pad_summary)

TREAT_LABEL_COLS <- c(
  "Non-Polluted + Ambient" = unname(TREAT_FILL["NP_Ambient"]),
  "Non-Polluted + Warm"    = unname(TREAT_FILL["NP_Warmed"]),
  "Polluted + Ambient"     = unname(TREAT_FILL["P_Ambient"]),
  "Polluted + Warm"        = unname(TREAT_FILL["P_Warmed"])
)

raw_plot <- ggplot(pad_summary, aes(x = Treatment, y = AvgAlgae, fill = Treatment)) +
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
  scale_y_continuous(
    name = "% macroalgal cover",
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    title = "Macroalgal cover by Warming × Pollution treatment",
    x = NULL,
    fill = "Treatment"
  ) +
  theme_minimal() +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(raw_plot)

## 2. Model variants & diagnostics --------------------------------------------
#   Summarise evidence for the modelling decisions used in main code.
## 2.1 Final model + diagnostics
algae_beta_final <- glmmTMB(
  prop_algae ~ Warming * Pollution +
    (1 + Warming * Pollution | `0-Country/site`),
  data   = pad_summary,
  family = beta_family(link = "logit"),
  control = glmmTMBControl(
    optimizer = optim,
    optArgs   = list(method = "BFGS"),
    profile   = TRUE
  )
)

print(summary(algae_beta_final))
print(make_diag_panel_glmmTMB(algae_beta_final, pad_summary, "Final beta–logit"))

## 2.2 Alternative variants
# No optimiser
models <- list(
  `Beta–logit` = glmmTMB(
    prop_algae ~ Warming * Pollution + (1 + Warming * Pollution | `0-Country/site`),
    data = pad_summary, family = beta_family(link = "logit") # Convergence issues
  ),
  # cloglog
  `Beta–cloglog` = glmmTMB(
    prop_algae ~ Warming * Pollution + (1 + Warming * Pollution | `0-Country/site`),
    data = pad_summary, family = beta_family(link = "cloglog") # Convergence issues
  ),
  # zero inflation
  `Zero-inflated beta` = glmmTMB(
    prop_algae ~ Warming * Pollution + (1 + Warming * Pollution | `0-Country/site`),
    ziformula = ~1,
    data = pad_summary, family = beta_family(link = "logit")
  ) # Convergence issues
)

for (nm in names(models)) {
  cat("Variant:", nm, "\n")
  print(summary(models[[nm]]))
  print(make_diag_panel_glmmTMB(models[[nm]], pad_summary, nm))
}

## 3. Site-level insights: intercept BLUPs --------------------------
re_int   <- ranef(algae_beta_final, condVar = TRUE)$cond$`0-Country/site`
pv_array <- attr(re_int, "condVar")
se_int   <- safe_sqrt(pv_array[1, 1, ])

site_blups <- tibble(
  Site       = clean_site(rownames(re_int)),
  blup_logit = re_int[, "(Intercept)"],
  se_logit   = se_int
)

global_logit <- unname(fixef(algae_beta_final)$cond["(Intercept)"])
global_prop  <- plogis(global_logit)

# pct_dev here is percentage-point deviation in cover
site_pct <- site_blups %>%
  mutate(
    logit_site  = global_logit + blup_logit,
    p_baseline  = plogis(logit_site),
    pct_dev     = (p_baseline - global_prop) * 100,
    lower_logit = logit_site - 1.96 * se_logit,
    upper_logit = logit_site + 1.96 * se_logit,
    lower_prop  = plogis(lower_logit),
    upper_prop  = plogis(upper_logit),
    lower_pct   = (lower_prop - global_prop) * 100,
    upper_pct   = (upper_prop - global_prop) * 100
  ) %>%
  left_join(
    pad_summary %>%
      distinct(`0-Country/site`, .keep_all = TRUE) %>%
      transmute(
        Site         = clean_site(`0-Country/site`),
        Latitude     = as.numeric(`0-Latitude`),
        Is_temperate = `4-Is_temperate`
      ),
    by = "Site"
  )

cat_plot <- ggplot(site_pct, aes(x = pct_dev, y = reorder(Site, pct_dev))) +
  geom_errorbarh(aes(xmin = lower_pct, xmax = upper_pct),
                 height = 0.2, color = "gray60") +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Site-level baseline deviations (beta–logit model)",
    subtitle = "Percentage-point deviation in baseline % cover from global mean",
    x        = "Deviation from global mean cover (percentage-points)",
    y        = NULL
  ) +
  theme_minimal()

print(cat_plot)

## 4. Site-level stressor effects (caterpillar plot on % odds scale; main code shows on model scale) --------------------

re_site  <- ranef(algae_beta_final, condVar = TRUE)$cond$`0-Country/site`
cv_array <- attr(re_site, "condVar")

site_slopes <- tibble(
  Site_raw         = rownames(re_site),
  blup_Warming     = re_site[, "Warming"],
  blup_Pollution   = re_site[, "Pollution"],
  blup_Interaction = re_site[, "Warming:Pollution"],
  se_Warming       = safe_sqrt(cv_array[2, 2, ]),
  se_Pollution     = safe_sqrt(cv_array[3, 3, ]),
  se_Interaction   = safe_sqrt(cv_array[4, 4, ])
) %>%
  mutate(Site = clean_site(Site_raw)) %>%
  pivot_longer(
    cols         = starts_with("blup_"),
    names_to     = "Term",
    names_prefix = "blup_",
    values_to    = "blup_link"
  ) %>%
  mutate(
    se_link = case_when(
      Term == "Warming"     ~ se_Warming,
      Term == "Pollution"   ~ se_Pollution,
      Term == "Interaction" ~ se_Interaction
    ),
    Term = factor(Term, levels = c("Warming", "Pollution", "Interaction"))
  ) %>%
  select(Site, Term, blup_link, se_link)

fe <- fixef(algae_beta_final)$cond
vc <- vcov(algae_beta_final)$cond

all_df <- site_slopes %>%
  mutate(
    beta_global = case_when(
      Term == "Warming"     ~ fe["Warming"],
      Term == "Pollution"   ~ fe["Pollution"],
      Term == "Interaction" ~ fe["Warming:Pollution"]
    ),
    slope_link = beta_global + blup_link,
    lower_link = slope_link - 1.96 * se_link,
    upper_link = slope_link + 1.96 * se_link,
    pct_odds   = (exp(slope_link) - 1) * 100,
    lower_pct  = (exp(lower_link) - 1) * 100,
    upper_pct  = (exp(upper_link) - 1) * 100
  )

global_df <- tibble(
  Term       = factor(c("Warming", "Pollution", "Interaction"),
                      levels = c("Warming", "Pollution", "Interaction")),
  slope_link = c(fe["Warming"], fe["Pollution"], fe["Warming:Pollution"]),
  se_link    = c(
    safe_sqrt(vc["Warming", "Warming"]),
    safe_sqrt(vc["Pollution", "Pollution"]),
    safe_sqrt(vc["Warming:Pollution", "Warming:Pollution"])
  )
) %>%
  mutate(
    lower_link = slope_link - 1.96 * se_link,
    upper_link = slope_link + 1.96 * se_link,
    pct_odds   = (exp(slope_link) - 1) * 100,
    lower_pct  = (exp(lower_link) - 1) * 100,
    upper_pct  = (exp(upper_link) - 1) * 100,
    Site       = "<b>Global fixed effect</b>"
  )

site_info <- pad_summary %>%
  distinct(`0-Country/site`, .keep_all = TRUE) %>%
  transmute(
    Site     = clean_site(`0-Country/site`),
    Latitude = as.numeric(`0-Latitude`)
  )

ord_sites <- site_info %>%
  semi_join(all_df %>% distinct(Site), by = "Site") %>%
  arrange(desc(Latitude)) %>%
  pull(Site)

y_levels <- rev(c("<b>Global fixed effect</b>", ord_sites))

p_pct_alg <- ggplot() +
  geom_rect(
    data  = global_df,
    aes(xmin = lower_pct, xmax = upper_pct),
    ymin  = -Inf, ymax = Inf,
    fill  = "grey80", alpha = 0.3
  ) +
  geom_vline(
    data = global_df,
    aes(xintercept = pct_odds),
    linetype = "dotted",
    colour   = "black"
  ) +
  geom_point(
    data   = global_df,
    aes(x = pct_odds, y = Site, fill = Term),
    shape  = 23,
    size   = 5,
    colour = "black"
  ) +
  geom_errorbarh(
    data   = all_df,
    aes(y = Site, xmin = lower_pct, xmax = upper_pct),
    height = 0.2,
    colour = "#666666"
  ) +
  geom_point(
    data = all_df,
    aes(x = pct_odds, y = Site, colour = Term),
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
  scale_x_continuous(
    name   = "% change in odds of algal cover",
    labels = function(x) paste0(x, "%")
  ) +
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

print(p_pct_alg)