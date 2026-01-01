###############################################################################
# Structural Equation Model (main analyses)
#   Data prep to pad-level SEM dataset
#   Final reduced / harmonised submodels (no W×P interactions)
#   Primary SEM + sensitivity SEM including Barnacles %~~% Grazers
#   Reportables: Fisher’s C, d-sep claims, standardised coefs, R²
#   Effects tables (direct / indirect / total), with optional bootstrap CIs
#
# *Additional code contains:
#   A priori interaction SEM and its diagnostics
#   Grazer model sensitivity/outlier checks
#   Full vs reduced SEM comparisons using AIC
#   Common-site sensitivity (intersection across submodels)
#   Random-effects structure alternatives (RI-only)
#   Competition-as-covariance alternative (Algae %~~% Barnacles)
###############################################################################

## 0. Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(DHARMa)
library(gt)
library(scales)
library(glue)
library(tibble)

set.seed(1)

## 0.1 Global constants & toggles

# NOTE: beta regression is not supported directly by piecewiseSEM.
# We therefore clip algae proportions away from 0/1 and use logit(p) as a
# pragmatic approximation of the beta–logit GLMM link.
CLIP_LO <- 0.001
CLIP_HI <- 0.999

# Bootstrap settings for indirect effect CIs:
#   nboot = 1000 is the setting used (~20-30mins run time).
#   For quick functionality checks, set NBOOT_TEST where indicated for shorter processing time.
RUN_BOOTSTRAP <- TRUE # Alternatively, if just checking code works and no need for bootstrapping, set this 'RUN_BOOTSTRAP' to 'FALSE'
NBOOT_FINAL   <- 1000
NBOOT_TEST    <- 100

## 0.2 Helper functions

# Keep proportions strictly within (0,1) before logit transform
clip01 <- function(p, lo = CLIP_LO, hi = CLIP_HI) pmin(pmax(p, lo), hi)

# Column selectors for piecewiseSEM::coefs() output
.pick_col <- function(df, patt, not = NULL) {
  ix <- grep(patt, names(df), ignore.case = TRUE)
  if (!is.null(not)) ix <- setdiff(ix, grep(not, names(df), ignore.case = TRUE))
  if (!length(ix)) stop("Couldn't find column matching: ", patt)
  names(df)[ix[1]]
}

# Normalise variable labels returned by coefs() (drops trailing digits)
.norm <- function(x) gsub("\\d+$", "", x)

# Extract standardised effect for a given (response, predictor)
.get_std <- function(sem_obj, resp, pred) {
  cf <- as_tibble(piecewiseSEM::coefs(sem_obj, standardize = "scale"),
                  .name_repair = "unique")
  resp_col <- .pick_col(cf, "^Resp")
  pred_col <- .pick_col(cf, "^Pred")
  std_col  <- .pick_col(cf, "Std(\\.|_)?Est|^std", not = "error|se")
  
  hit <- cf %>%
    mutate(
      Resp = .norm(.data[[resp_col]]),
      Pred = .norm(.data[[pred_col]])
    ) %>%
    filter(Resp == resp, Pred == pred) %>%
    pull(std_col)
  
  if (!length(hit)) return(0)
  suppressWarnings(as.numeric(hit[1]))
}

# Extract raw (unstandardised link-scale) effect for a given (response, predictor)
.get_raw <- function(sem_obj, resp, pred) {
  cf <- as_tibble(piecewiseSEM::coefs(sem_obj, standardize = "none"),
                  .name_repair = "unique")
  resp_col <- .pick_col(cf, "^Resp")
  pred_col <- .pick_col(cf, "^Pred")
  
  # piecewiseSEM can sometimes name the unstandardised coefficient column differently, correct
  cand <- grep("^(Est(\\.|imate)?|coef|b)$|^Estimate$",
               names(cf), ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, grep("Std", names(cf), ignore.case = TRUE, value = TRUE))
  est_col <- if (length(cand)) cand[1] else "Estimate"
  if (!est_col %in% names(cf)) return(0)
  
  hit <- cf %>%
    mutate(
      Resp = .norm(.data[[resp_col]]),
      Pred = .norm(.data[[pred_col]])
    ) %>%
    filter(Resp == resp, Pred == pred) %>%
    pull(est_col)
  
  if (!length(hit)) return(0)
  suppressWarnings(as.numeric(hit[1]))
}

# Build direct/indirect/total effects table (standardised)
effects_table_std <- function(sem_obj) {
  W <- "Warming"; P <- "Pollution"
  A <- "logit_algae"; B <- "log_AvgBarnacles"; G <- "log_AvgGrazers"
  
  # Direct links (standardised)
  d_WA <- .get_std(sem_obj, A, W)
  d_WB <- .get_std(sem_obj, B, W)
  d_WG <- 0  # no direct Warming -> Grazers path in the DAG
  
  d_PA <- .get_std(sem_obj, A, P)
  d_PB <- .get_std(sem_obj, B, P)
  d_PG <- .get_std(sem_obj, G, P)
  
  d_AB <- .get_std(sem_obj, B, A)
  d_AG <- .get_std(sem_obj, G, A)
  
  # Indirect via algae (product of link coefficients)
  i_WB <- d_WA * d_AB
  i_WG <- d_WA * d_AG
  i_PB <- d_PA * d_AB
  i_PG <- d_PA * d_AG
  
  tibble(
    From = c(W,  W,  W,  P,  P,  P,  A,  A),
    To   = c(A,  B,  G,  A,  B,  G,  B,  G),
    Type = "Direct",
    Estimate = c(d_WA, d_WB, d_WG, d_PA, d_PB, d_PG, d_AB, d_AG)
  ) %>%
    bind_rows(
      tibble(
        From = c(W,  W,  P,  P),
        To   = c(B,  G,  B,  G),
        Type = "Indirect (via algae)",
        Estimate = c(i_WB, i_WG, i_PB, i_PG),
        Indirect_Path = c(
          "Warming → Algae → Barnacles",
          "Warming → Algae → Grazers",
          "Pollution → Algae → Barnacles",
          "Pollution → Algae → Grazers"
        )
      )
    ) %>%
    mutate(
      Total = dplyr::case_when(
        Type == "Direct" ~ Estimate,
        Type != "Direct" & From == W & To == B ~ i_WB + d_WB,
        Type != "Direct" & From == W & To == G ~ i_WG + d_WG,
        Type != "Direct" & From == P & To == B ~ i_PB + d_PB,
        Type != "Direct" & From == P & To == G ~ i_PG + d_PG,
        TRUE ~ Estimate
      ),
      Type = factor(Type, levels = c("Direct", "Indirect (via algae)"))
    ) %>%
    arrange(To, Type, From)
}

# Build direct/indirect/total effects table (raw link-scale)
effects_table_raw <- function(sem_obj) {
  W <- "Warming"; P <- "Pollution"
  A <- "logit_algae"; B <- "log_AvgBarnacles"; G <- "log_AvgGrazers"
  
  d_WA <- .get_raw(sem_obj, A, W)
  d_WB <- .get_raw(sem_obj, B, W)
  d_WG <- 0
  
  d_PA <- .get_raw(sem_obj, A, P)
  d_PB <- .get_raw(sem_obj, B, P)
  d_PG <- .get_raw(sem_obj, G, P)
  
  d_AB <- .get_raw(sem_obj, B, A)
  d_AG <- .get_raw(sem_obj, G, A)
  
  i_WB <- d_WA * d_AB
  i_WG <- d_WA * d_AG
  i_PB <- d_PA * d_AB
  i_PG <- d_PA * d_AG
  
  tibble(
    From = c(W,  W,  W,  P,  P,  P,  A,  A),
    To   = c(A,  B,  G,  A,  B,  G,  B,  G),
    Type = "Direct",
    Estimate = c(d_WA, d_WB, d_WG, d_PA, d_PB, d_PG, d_AB, d_AG)
  ) %>%
    bind_rows(
      tibble(
        From = c(W,  W,  P,  P),
        To   = c(B,  G,  B,  G),
        Type = "Indirect (via algae)",
        Estimate = c(i_WB, i_WG, i_PB, i_PG),
        Indirect_Path = c(
          "Warming → Algae → Barnacles",
          "Warming → Algae → Grazers",
          "Pollution → Algae → Barnacles",
          "Pollution → Algae → Grazers"
        )
      )
    ) %>%
    mutate(
      Total = dplyr::case_when(
        Type == "Direct" ~ Estimate,
        Type != "Direct" & From == W & To == B ~ i_WB + d_WB,
        Type != "Direct" & From == W & To == G ~ i_WG + d_WG,
        Type != "Direct" & From == P & To == B ~ i_PB + d_PB,
        Type != "Direct" & From == P & To == G ~ i_PG + d_PG,
        TRUE ~ Estimate
      ),
      Type = factor(Type, levels = c("Direct", "Indirect (via algae)"))
    ) %>%
    arrange(To, Type, From)
}

# Extract direct-path p-values from fitted SEM (for table annotation)
get_direct_pvals <- function(sem_obj) {
  co_raw <- piecewiseSEM::coefs(sem_obj, standardize = "none")
  co <- tibble::as_tibble(co_raw, .name_repair = "unique")
  
  pred_col <- .pick_col(co, "^Pred")
  resp_col <- .pick_col(co, "^Resp")
  p_col    <- .pick_col(co, "P(\\.|_|$)|P.value")
  
  co %>%
    transmute(
      From = .norm(.data[[pred_col]]),
      To   = .norm(.data[[resp_col]]),
      p    = suppressWarnings(as.numeric(.data[[p_col]]))
    ) %>%
    filter(!is.na(From), !is.na(To), From != "(Intercept)") %>%
    distinct(From, To, .keep_all = TRUE)
}

# Table formatting helpers
to_num <- function(x) suppressWarnings(as.numeric(gsub("\\s|\\+", "", as.character(x))))

ci_excludes_zero <- function(lo, hi) {
  lo_n <- to_num(lo)
  hi_n <- to_num(hi)
  !is.na(lo_n) & !is.na(hi_n) & (lo_n * hi_n > 0)
}

nice_ci <- function(est, lo, hi, digits = 3) {
  sprintf(
    "%+.*f (%+.*f, %+.*f)",
    digits, to_num(est),
    digits, to_num(lo),
    digits, to_num(hi)
  )
}

flag_mark <- function(x) ifelse(is.na(x), "—", ifelse(as.logical(x), "✓", "✗"))

build_sem_table <- function(effects_summary, direct_p, title = NULL, subtitle = NULL) {
  
  notes_text <- glue(
    "**Notes.** Direct paths use the mixed-model p-values (`p (direct)`). ",
    "Bootstrap CIs are percentile 95% CIs from the parametric SEM bootstrap. ",
    "Markers: Y = criterion met; N = criterion not met; NA = not applicable."
  )
  
  effects_clean <- effects_summary %>%
    mutate(
      Effect   = ifelse(grepl("Indirect", Type, ignore.case = TRUE),
                        "Indirect (via algae)", "Direct"),
      Estimate = to_num(Estimate),
      Est_low  = to_num(Est_low),
      Est_high = to_num(Est_high),
      Total    = to_num(Total),
      Tot_low  = to_num(Tot_low),
      Tot_high = to_num(Tot_high)
    ) %>%
    left_join(direct_p, by = c("From", "To")) %>%
    mutate(
      Effect_CI = nice_ci(Estimate, Est_low, Est_high),
      Total_CI  = nice_ci(Total,    Tot_low, Tot_high),
      p_direct  = ifelse(
        Effect == "Direct" & !is.na(p),
        scales::pvalue(p, accuracy = 0.001),
        "—"
      ),
      p_sig     = Effect == "Direct" & !is.na(p) & p < 0.05,
      ci_est    = ci_excludes_zero(Est_low, Est_high),
      ci_tot    = ci_excludes_zero(Tot_low, Tot_high)
    ) %>%
    transmute(
      To, From, Effect,
      Path = dplyr::if_else(
        Effect == "Indirect (via algae)",
        Indirect_Path,
        NA_character_
      ),
      Effect_CI, Total_CI, p_direct,
      p_sig, ci_est, ci_tot
    ) %>%
    arrange(To, factor(Effect, levels = c("Direct", "Indirect (via algae)")), From)
  
  tbl <- effects_clean %>%
    mutate(
      p_sig_flag  = flag_mark(p_sig),
      ci_est_flag = flag_mark(ci_est),
      ci_tot_flag = flag_mark(ci_tot)
    ) %>%
    select(
      To, From, Effect, Path,
      Effect_CI, Total_CI, p_direct,
      p_sig_flag, ci_est_flag, ci_tot_flag
    ) %>%
    gt(rowname_col = NULL, groupname_col = "To") %>%
    cols_label(
      From        = "From",
      Effect      = "Type",
      Path        = "Path",
      Effect_CI   = "Estimate (95% CI)",
      Total_CI    = "Total (95% CI)",
      p_direct    = "p (direct)",
      p_sig_flag  = "p < 0.05",
      ci_est_flag = "CI≠0 (effect)",
      ci_tot_flag = "CI≠0 (total)"
    ) %>%
    tab_spanner(
      label = "Evidence of effect",
      columns = c(p_direct, p_sig_flag, ci_est_flag, ci_tot_flag)
    ) %>%
    cols_align(
      align = "center",
      columns = c(p_direct, p_sig_flag, ci_est_flag, ci_tot_flag)
    ) %>%
    fmt_missing(everything(), missing_text = "—") %>%
    tab_options(
      table.font.names = "Helvetica",
      table.font.size  = 12
    ) %>%
    opt_row_striping() %>%
    tab_source_note(md(notes_text))
  
  if (!is.null(title) || !is.null(subtitle)) {
    tbl <- tab_header(
      tbl,
      title    = md(title %||% ""),
      subtitle = md(subtitle %||% "")
    )
  }
  
  tbl
}

## 0.3 Parametric bootstrap helpers
#
#   Simulate a new response vector from each fitted submodel (algae, barnacles, grazers)
#   Refit each submodel to its simulated response (keeping the same RE structure)
#   Rebuild SEM (including Barnacles %~~% Grazers covariance)
#   Recompute direct/indirect/total effects
#   Summarise percentile CIs across replicates
#
#   This is computationally heavy. The script prints elapsed time for chosen nboot; 
#   For quick checks, use small nboot first (e.g. NBOOT_TEST).

.resp_name <- function(mod) {
  f <- formula(mod)
  as.character(f[[2]])
}

.rebuild_with_y <- function(mod, y) {
  d <- model.frame(mod)
  d[[.resp_name(mod)]] <- as.numeric(y)
  
  ref <- try(
    suppressWarnings(update(mod, data = d, control = lmerControl(calc.derivs = FALSE))),
    silent = TRUE
  )
  if (inherits(ref, "try-error")) {
    ref <- try(suppressWarnings(update(mod, data = d)), silent = TRUE)
  }
  if (inherits(ref, "try-error")) return(NULL)
  ref
}

boot_sem_both <- function(
    nboot = 1000,
    seed = 1,
    algae_mod,
    barn_mod,
    graz_mod
) {
  set.seed(seed)
  
  # Prototypes keep row order stable even if a replicate fails to return some rows
  proto_std <- effects_table_std(sem_BGcorr) %>%
    transmute(From, To, Type,
              Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path))
  
  proto_raw <- effects_table_raw(sem_BGcorr) %>%
    transmute(From, To, Type,
              Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path))
  
  do_one <- function(i) {
    
    # Simulate responses
    yA <- simulate(algae_mod, nsim = 1)[[1]]
    yB <- simulate(barn_mod,  nsim = 1)[[1]]
    yG <- simulate(graz_mod,  nsim = 1)[[1]]
    
    # Refit submodels on simulated responses
    refA <- .rebuild_with_y(algae_mod, yA); if (is.null(refA)) return(NULL)
    refB <- .rebuild_with_y(barn_mod,  yB); if (is.null(refB)) return(NULL)
    refG <- .rebuild_with_y(graz_mod,  yG); if (is.null(refG)) return(NULL)
    
    # Rebuild SEM (including residual covariance used in the paper)
    sem_i <- try(
      piecewiseSEM::psem(
        refA, refB, refG,
        log_AvgBarnacles %~~% log_AvgGrazers
      ),
      silent = TRUE
    )
    if (inherits(sem_i, "try-error")) return(NULL)
    
    # Extract effects (std + raw), and right_join prototypes for stability
    eff_std <- effects_table_std(sem_i) %>%
      mutate(Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path)) %>%
      select(From, To, Type, Indirect_Path, Estimate, Total) %>%
      right_join(proto_std, by = c("From", "To", "Type", "Indirect_Path")) %>%
      arrange(To, Type, From) %>%
      mutate(Indirect_Path = ifelse(Indirect_Path == "__NA__", NA_character_, Indirect_Path))
    
    eff_raw <- effects_table_raw(sem_i) %>%
      mutate(Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path)) %>%
      select(From, To, Type, Indirect_Path, Estimate, Total) %>%
      right_join(proto_raw, by = c("From", "To", "Type", "Indirect_Path")) %>%
      arrange(To, Type, From) %>%
      mutate(Indirect_Path = ifelse(Indirect_Path == "__NA__", NA_character_, Indirect_Path))
    
    list(std = eff_std, raw = eff_raw)
  }
  
  pb <- txtProgressBar(min = 0, max = nboot, style = 3)
  res <- vector("list", nboot)
  for (i in seq_len(nboot)) {
    res[[i]] <- do_one(i)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ok <- purrr::compact(res)
  if (!length(ok)) stop("No successful bootstrap replicates.")
  
  boot_df_std <- bind_rows(lapply(ok, `[[`, "std"), .id = "replicate")
  boot_df_raw <- bind_rows(lapply(ok, `[[`, "raw"), .id = "replicate")
  
  sum_std <- boot_df_std %>%
    group_by(From, To, Type, Indirect_Path) %>%
    summarise(
      n        = n(),
      Est_low  = quantile(Estimate, probs = 0.025, na.rm = TRUE),
      Est_high = quantile(Estimate, probs = 0.975, na.rm = TRUE),
      Tot_low  = quantile(Total,    probs = 0.025, na.rm = TRUE),
      Tot_high = quantile(Total,    probs = 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
  
  sum_raw <- boot_df_raw %>%
    group_by(From, To, Type, Indirect_Path) %>%
    summarise(
      n        = n(),
      Est_low  = quantile(Estimate, probs = 0.025, na.rm = TRUE),
      Est_high = quantile(Estimate, probs = 0.975, na.rm = TRUE),
      Tot_low  = quantile(Total,    probs = 0.025, na.rm = TRUE),
      Tot_high = quantile(Total,    probs = 0.975, na.rm = TRUE),
      .groups  = "drop"
    )
  
  list(
    boot_draws_std   = boot_df_std,
    boot_summary_std = sum_std,
    boot_draws_raw   = boot_df_raw,
    boot_summary_raw = sum_raw
  )
}

## 1. Data preparation & pad-level SEM dataset --------------------------------

# Load raw long data
hotmess_full <- readr::read_csv(
  here::here("Dataframes", "Full global HotMess data.csv"),
  show_col_types = FALSE
)

# Treatments as 0/1 (control = 0)
hotmess_full <- hotmess_full %>%
  mutate(
    Warming   = if_else(`0-Pad colour` == "B", 1L, 0L),
    Pollution = if_else(`0-Pad region` == "P", 1L, 0L)
  )

# Combine grazer counts (surface + sides), as with GLMM
hotmess_full <- hotmess_full %>%
  mutate(
    `3-Grazer total count` =
      `3-Grazer total count on pad surface` +
      `3-Grazer total count on pad sides`
  ) %>%
  dplyr::select(
    -`3-Grazer total count on pad surface`,
    -`3-Grazer total count on pad sides`
  )

# Summarise to pad level (average across timepoints)
sem_data <- hotmess_full %>%
  group_by(`0-Country/site`, `0-Pad ID`) %>%
  summarise(
    `0-Pad region`   = first(`0-Pad region`),
    `0-Pad colour`   = first(`0-Pad colour`),
    `0-Latitude`     = first(`0-Latitude`),
    `4-Is_temperate` = first(`4-Is_temperate`),
    Warming          = first(Warming),
    Pollution        = first(Pollution),
    
    AvgBarnacles = mean(`3-Barnacle total count`,       na.rm = TRUE),
    AvgAlgae     = mean(`3-Macroalgae total pct cover`, na.rm = TRUE),
    AvgGrazers   = mean(`3-Grazer total count`,         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(AvgBarnacles))

# Transformations used by the SEM:
#   algae: logit(clip(prop))
#   barnacle + grazer abundance: log(count + 1)
sem_data <- sem_data %>%
  mutate(
    site             = factor(`0-Country/site`),
    prop_algae       = AvgAlgae / 100,
    logit_algae      = qlogis(clip01(prop_algae)),
    log_AvgBarnacles = log(AvgBarnacles + 1),
    log_AvgGrazers   = log(AvgGrazers + 1)
  )

# Per-submodel site filtering:
#   algae: drop Chile (macroalgae <1%)
#   barnacles: drop Argentina (0 barnacles)
#   grazers: keep all sites
dat_algae     <- sem_data %>% filter(site != "Chile")     %>% droplevels()
dat_barnacles <- sem_data %>% filter(site != "Argentina") %>% droplevels()
dat_grazers   <- sem_data                                  %>% droplevels()

## 2. Final submodels for SEM (reduced, harmonised) ----------------------------

# Model design (final, paper-facing):
#   Algae:     logit_algae ~ Warming + Pollution
#   Barnacles: log_AvgBarnacles ~ Warming + Pollution + logit_algae
#   Grazers:   log_AvgGrazers   ~ Pollution + logit_algae
#
# Random effects:
#   (1 + Pollution | site) for all three equations
#   Represents site heterogeneity in the pollution effect.

algae_mod_red <- lmer(
  logit_algae ~ Warming + Pollution + (1 + Pollution | site),
  data = dat_algae
)

barn_mod_red <- lmer(
  log_AvgBarnacles ~ Warming + Pollution + logit_algae + (1 + Pollution | site),
  data = dat_barnacles
)

graz_mod_logit <- lmer(
  log_AvgGrazers ~ logit_algae + Pollution + (1 + Pollution | site),
  data = dat_grazers
)

# Quick singularity verification
print(isSingular(algae_mod_red, tol = 1e-5))
print(isSingular(barn_mod_red,  tol = 1e-5))
print(isSingular(graz_mod_logit, tol = 1e-5))

## 3. SEM construction ------------------------------------------

# 3.1 Primary reduced SEM (no residual covariance between barnacles and grazers)
sem_red <- psem(algae_mod_red, barn_mod_red, graz_mod_logit)

print(summary(sem_red))
fc_red <- fisherC(sem_red)
print(fc_red)

# Fisher’s C p-value
C_red  <- unname(as.numeric(fc_red[[1]]))
df_red <- unname(as.numeric(fc_red[[2]]))
P_red  <- pchisq(C_red, df_red, lower.tail = FALSE)

print(data.frame(Model = "Primary reduced", FishersC = C_red, df = df_red, Pvalue = P_red))
print(coefs(sem_red, standardize = "scale"))
print(rsquared(sem_red))

# 3.2 Sensitivity SEM: Barnacles %~~% Grazers
#
# Interpretation:
#   This adds a residual covariance (not an a priori directed path) to absorb
#   an unmodelled shared driver (e.g. microhabitat, local food, etc.).

sem_BGcorr <- psem(
  algae_mod_red,
  barn_mod_red,
  graz_mod_logit,
  log_AvgBarnacles %~~% log_AvgGrazers
)
# Note, warning refers to the LL-based Chi-squared (which is not the primary fit statistic we report).
# Fisher’s C + d-sep claims remain the appropriate global fit outputs here.

print(summary(sem_BGcorr))
fc_BG <- fisherC(sem_BGcorr)
print(fc_BG)

C_bg  <- unname(as.numeric(fc_BG[[1]]))
df_bg <- unname(as.numeric(fc_BG[[2]]))
P_bg  <- pchisq(C_bg, df_bg, lower.tail = FALSE)

print(data.frame(Model = "Reduced + Barnacles%~~%Grazers", FishersC = C_bg, df = df_bg, Pvalue = P_bg))
print(coefs(sem_BGcorr, standardize = "scale"))
print(rsquared(sem_BGcorr))

## 4. Effects tables (point estimates + bootstrap CIs)

# We treat sem_BGcorr as the reporting SEM here
effects_point_std <- effects_table_std(sem_BGcorr) %>%
  mutate(Estimate = round(Estimate, 3), Total = round(Total, 3))

effects_point_raw <- effects_table_raw(sem_BGcorr) %>%
  mutate(Estimate = round(Estimate, 3), Total = round(Total, 3))

# Bootstrap CIs
if (isTRUE(RUN_BOOTSTRAP)) {
  
  nboot_use <- NBOOT_FINAL
  # If just checking that the bootstrap machinery still works after edits,
  # switch to NBOOT_TEST by changing the line below.
  # nboot_use <- NBOOT_TEST
  
  cat("\n====================\nBOOTSTRAP\n====================\n")
  cat("Running parametric SEM bootstrap with nboot =", nboot_use, "\n")
  cat("This prints elapsed time when complete. For quick checks, reduce nboot.\n")
  
  set.seed(42)
  t_boot <- system.time({
    boot_out <- boot_sem_both(
      nboot = nboot_use,
      seed  = 2025,
      algae_mod = algae_mod_red,
      barn_mod  = barn_mod_red,
      graz_mod  = graz_mod_logit
    )
  })
  
  cat("\nBootstrap elapsed time:\n")
  print(t_boot)
  
  # Combine point estimates + bootstrap CIs
  effects_summary_std <- effects_point_std %>%
    left_join(
      boot_out$boot_summary_std,
      by = c("From", "To", "Type", "Indirect_Path")
    ) %>%
    mutate(
      Est_low  = round(Est_low,  3),
      Est_high = round(Est_high, 3),
      Tot_low  = round(Tot_low,  3),
      Tot_high = round(Tot_high, 3)
    ) %>%
    select(
      From, To, Type, Indirect_Path,
      Estimate, Est_low, Est_high,
      Total,    Tot_low, Tot_high
    )
  
  effects_summary_raw <- effects_point_raw %>%
    left_join(
      boot_out$boot_summary_raw,
      by = c("From", "To", "Type", "Indirect_Path")
    ) %>%
    mutate(
      Est_low  = round(Est_low,  3),
      Est_high = round(Est_high, 3),
      Tot_low  = round(Tot_low,  3),
      Tot_high = round(Tot_high, 3)
    ) %>%
    select(
      From, To, Type, Indirect_Path,
      Estimate, Est_low, Est_high,
      Total,    Tot_low, Tot_high
    )
  
} else {
  
  # If bootstrap is OFF, we still return tables with placeholder CI columns
  effects_summary_std <- effects_point_std %>%
    mutate(Est_low = NA_real_, Est_high = NA_real_, Tot_low = NA_real_, Tot_high = NA_real_) %>%
    select(From, To, Type, Indirect_Path, Estimate, Est_low, Est_high, Total, Tot_low, Tot_high)
  
  effects_summary_raw <- effects_point_raw %>%
    mutate(Est_low = NA_real_, Est_high = NA_real_, Tot_low = NA_real_, Tot_high = NA_real_) %>%
    select(From, To, Type, Indirect_Path, Estimate, Est_low, Est_high, Total, Tot_low, Tot_high)
}

# Build tables (standardised + raw)
direct_p <- get_direct_pvals(sem_BGcorr)

tbl_gt_std <- build_sem_table(
  effects_summary = effects_summary_std,
  direct_p        = direct_p,
  title    = "**SEM effects — Standardised coefficients**",
  subtitle = "Direct + indirect (via algae) + total effects; CIs are bootstrap percentile CIs if RUN_BOOTSTRAP=TRUE."
)

tbl_gt_raw <- build_sem_table(
  effects_summary = effects_summary_raw,
  direct_p        = direct_p,
  title    = "**SEM effects — Raw (link-scale) coefficients**",
  subtitle = "Raw link-scale effects; indirect = product of raw links; CIs are bootstrap percentile CIs if RUN_BOOTSTRAP=TRUE."
)

tbl_gt_std
tbl_gt_raw