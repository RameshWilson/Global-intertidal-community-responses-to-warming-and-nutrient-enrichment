###############################################################################
# Structural Equation Model (additional analyses)
#   A priori interaction SEM and diagnostics
#   Algae as Gaussian proportion (early attempt)
#   Grazer sensitivity checks (outliers, algae scale, RE structure)
#   Full vs reduced comparisons (including AIC; documented as “evidence” only)
#   Common-site sensitivity (intersection across submodels)
#   Alternative RE structure (RI-only)
#   Competition as residual covariance (Algae %~~% Barnacles)
#   DiagrammeR “visualisation
#   Bootstrap troubleshooting
###############################################################################

## 0. Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(DHARMa)
library(DiagrammeR)
library(glue)
library(scales)
library(tibble)
library(gt)

set.seed(1)

## 0.1 Global constants & toggles ---------------------------------------------

CLIP_LO <- 0.001
CLIP_HI <- 0.999

# DHARMa simulations per check
DHARMA_NSIM <- 1000

# For quick bootstrap sanity checks in ADDITIONAL
RUN_BOOT_SANITY <- FALSE
NBOOT_SANITY    <- 100

## 0.2 Helper functions

# Keep proportions strictly within (0,1) before logit transform
clip01 <- function(p, lo = CLIP_LO, hi = CLIP_HI) pmin(pmax(p, lo), hi)

# Diagnostics: DHARMa plot + key tests
dharma_quick <- function(mod, name = "Model", nsim = DHARMA_NSIM) {
  set.seed(1)
  sim <- DHARMa::simulateResiduals(mod, n = nsim)
  par(mfrow = c(1, 2)); plot(sim); par(mfrow = c(1, 1))
  cat("\n--- DHARMa tests:", name, "---\n")
  print(DHARMa::testUniformity(sim))
  print(DHARMa::testDispersion(sim))
  invisible(sim)
}

# Robust column pickers for piecewiseSEM::coefs() output
.pick_col <- function(df, patt, not = NULL) {
  ix <- grep(patt, names(df), ignore.case = TRUE)
  if (!is.null(not)) ix <- setdiff(ix, grep(not, names(df), ignore.case = TRUE))
  if (!length(ix)) stop("Couldn't find column matching: ", patt)
  names(df)[ix[1]]
}

.norm <- function(x) gsub("\\d+$", "", x)

.get_std <- function(sem_obj, resp, pred) {
  cf <- as_tibble(piecewiseSEM::coefs(sem_obj, standardize = "scale"),
                  .name_repair = "unique")
  resp_col <- .pick_col(cf, "^Resp")
  pred_col <- .pick_col(cf, "^Pred")
  std_col  <- .pick_col(cf, "Std(\\.|_)?Est|^std", not = "error|se")
  
  hit <- cf %>%
    mutate(Resp = .norm(.data[[resp_col]]),
           Pred = .norm(.data[[pred_col]])) %>%
    filter(Resp == resp, Pred == pred) %>%
    pull(std_col)
  
  if (!length(hit)) return(0)
  suppressWarnings(as.numeric(hit[1]))
}

.get_raw <- function(sem_obj, resp, pred) {
  cf <- as_tibble(piecewiseSEM::coefs(sem_obj, standardize = "none"),
                  .name_repair = "unique")
  resp_col <- .pick_col(cf, "^Resp")
  pred_col <- .pick_col(cf, "^Pred")
  
  cand <- grep("^(Est(\\.|imate)?|coef|b)$|^Estimate$",
               names(cf), ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, grep("Std", names(cf), ignore.case = TRUE, value = TRUE))
  est_col <- if (length(cand)) cand[1] else "Estimate"
  if (!est_col %in% names(cf)) return(0)
  
  hit <- cf %>%
    mutate(Resp = .norm(.data[[resp_col]]),
           Pred = .norm(.data[[pred_col]])) %>%
    filter(Resp == resp, Pred == pred) %>%
    pull(est_col)
  
  if (!length(hit)) return(0)
  suppressWarnings(as.numeric(hit[1]))
}

effects_table_std <- function(sem_obj) {
  W <- "Warming"; P <- "Pollution"
  A <- "logit_algae"; B <- "log_AvgBarnacles"; G <- "log_AvgGrazers"
  
  d_WA <- .get_std(sem_obj, A, W)
  d_WB <- .get_std(sem_obj, B, W)
  d_WG <- 0
  
  d_PA <- .get_std(sem_obj, A, P)
  d_PB <- .get_std(sem_obj, B, P)
  d_PG <- .get_std(sem_obj, G, P)
  
  d_AB <- .get_std(sem_obj, B, A)
  d_AG <- .get_std(sem_obj, G, A)
  
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

# Bootstrap helpers
.resp_name <- function(mod) as.character(formula(mod)[[2]])

.rebuild_with_y <- function(mod, y) {
  d <- model.frame(mod)
  d[[.resp_name(mod)]] <- as.numeric(y)
  ref <- try(suppressWarnings(update(mod, data = d, control = lmerControl(calc.derivs = FALSE))), silent = TRUE)
  if (inherits(ref, "try-error")) ref <- try(suppressWarnings(update(mod, data = d)), silent = TRUE)
  if (inherits(ref, "try-error")) return(NULL)
  ref
}

boot_sem_both <- function(nboot, seed, algae_mod, barn_mod, graz_mod, sem_template = NULL) {
  set.seed(seed)
  
  # Use template SEM object to define stable row order for returned tables
  if (is.null(sem_template)) stop("Provide sem_template so prototypes are well-defined.")
  
  proto_std <- effects_table_std(sem_template) %>%
    transmute(From, To, Type,
              Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path))
  
  proto_raw <- effects_table_raw(sem_template) %>%
    transmute(From, To, Type,
              Indirect_Path = ifelse(is.na(Indirect_Path), "__NA__", Indirect_Path))
  
  do_one <- function(i) {
    
    yA <- simulate(algae_mod, nsim = 1)[[1]]
    yB <- simulate(barn_mod,  nsim = 1)[[1]]
    yG <- simulate(graz_mod,  nsim = 1)[[1]]
    
    refA <- .rebuild_with_y(algae_mod, yA); if (is.null(refA)) return(NULL)
    refB <- .rebuild_with_y(barn_mod,  yB); if (is.null(refB)) return(NULL)
    refG <- .rebuild_with_y(graz_mod,  yG); if (is.null(refG)) return(NULL)
    
    sem_i <- try(
      piecewiseSEM::psem(refA, refB, refG, log_AvgBarnacles %~~% log_AvgGrazers),
      silent = TRUE
    )
    if (inherits(sem_i, "try-error")) return(NULL)
    
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
  
  list(std = boot_df_std, raw = boot_df_raw, n_ok = length(ok), n_try = nboot)
}

## 1. Data preparation & pad-level SEM dataset --------------------------------

# Load raw long data
hotmess_full <- readr::read_csv(
  here::here("Dataframes", "Full global HotMess data.csv"),
  show_col_types = FALSE
)

# Treatments as 0/1 (control = 0)
# Combine grazer counts (surface + sides), as with GLMM
hotmess_full <- hotmess_full %>%
  mutate(
    Warming   = if_else(`0-Pad colour` == "B", 1L, 0L),
    Pollution = if_else(`0-Pad region` == "P", 1L, 0L)
  ) %>%
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
  filter(!is.na(AvgBarnacles)) %>%
  mutate(
    site             = factor(`0-Country/site`),
    prop_algae       = AvgAlgae / 100,
    log_AvgBarnacles = log(AvgBarnacles + 1),
    log_AvgGrazers   = log(AvgGrazers + 1),
    logit_algae      = qlogis(clip01(prop_algae))
  )

# Per-submodel filtering
dat_algae     <- sem_data %>% filter(site != "Chile")     %>% droplevels()
dat_barnacles <- sem_data %>% filter(site != "Argentina") %>% droplevels()
dat_grazers   <- sem_data                                  %>% droplevels()

## 2. A priori submodels (with interactions; early attempt)

# A priori causal structure (early form):
#   Algae     ~ Warming * Pollution
#   Barnacles ~ Warming * Pollution + algae
#   Grazers   ~ algae + Pollution (no Warming effect likely; grazers can avoid hot pads)
#
# Original algae-on-proportion Gaussian approach; not used for results.
algae_mod_full_orig <- lmer(
  prop_algae ~ Warming * Pollution + (1 + Pollution | site),
  data = dat_algae
)

barn_mod_full_orig <- lmer(
  log_AvgBarnacles ~ Warming * Pollution + prop_algae + (1 + Pollution | site),
  data = dat_barnacles
)

graz_mod_full_orig <- lmer(
  log_AvgGrazers ~ prop_algae + Pollution + (1 + Pollution | site),
  data = dat_grazers
)

dharma_quick(algae_mod_full_orig, "Algae (prop, Gaussian; original)")
dharma_quick(barn_mod_full_orig,  "Barnacles (log+1; original)")
dharma_quick(graz_mod_full_orig,  "Grazers (log+1; original)")

## 3. Algae logit transform
# Motivation:
#   piecewiseSEM does not support beta regression, of which was employed for algae GLMMs
#   logit(clip(prop)) approximates a beta–logit link
#   diagnostics typically improve vs Gaussian on raw proportions

dat_algae <- dat_algae %>%
  mutate(
    prop_algae  = clip01(AvgAlgae / 100),
    logit_algae = qlogis(prop_algae)
  )

algae_mod_full <- lmer(
  logit_algae ~ Warming * Pollution + (1 + Pollution | site),
  data = dat_algae
)

dharma_quick(algae_mod_full, "Algae (logit clipped prop; interaction)")

## 4. Grazer sensitivity block (outliers / algae scale / slopes)

# We evidence that:
#   A few influential points can create mild non-uniformity;
#   Alternatives do not materially change inferences;
#   Therefore we keep the simpler a priori grazer structure in main

graz_mod_base <- lmer(
  log_AvgGrazers ~ prop_algae + Pollution + (1 + Pollution | site),
  data = dat_grazers
)

# 4.1 Sensitivity to top residual outliers (remove top 3 Pearson residuals)
prs  <- residuals(graz_mod_base, type = "pearson")
idx3 <- order(abs(prs), decreasing = TRUE)[1:3]
dat_grazers_sens <- dat_grazers[-idx3, ]

graz_mod_no3 <- lmer(
  log_AvgGrazers ~ prop_algae + Pollution + (1 + Pollution | site),
  data = dat_grazers_sens
)

print(summary(graz_mod_base)$coef)
print(summary(graz_mod_no3)$coef)

# 4.2 Alternative algae scale as predictor (logit instead of proportion)
dat_grazers <- dat_grazers %>%
  mutate(
    prop_algae_clip = clip01(prop_algae),
    logit_algae     = qlogis(prop_algae_clip)
  )

graz_mod_A <- lmer(
  log_AvgGrazers ~ logit_algae + Pollution + (1 + Pollution | site),
  data = dat_grazers
)

# 4.3 Add random slope for algae (uncorrelated slopes)
graz_mod_B <- lmer(
  log_AvgGrazers ~ prop_algae + Pollution +
    (1 | site) + (0 + Pollution | site) + (0 + prop_algae | site),
  data = dat_grazers
)

dharma_quick(graz_mod_base, "Grazers baseline")
dharma_quick(graz_mod_A,    "Grazers with logit_algae predictor")
dharma_quick(graz_mod_B,    "Grazers with extra random slope for algae")

print(AIC(graz_mod_base, graz_mod_A, graz_mod_B))
print(c(
  isSingular(graz_mod_base),
  isSingular(graz_mod_A),
  isSingular(graz_mod_B)
))

## 5. SEM construction: full vs reduced ------------------------

# 5.1 Full SEM (interactions + mixed algae scaling)

sem_full <- psem(
  algae_mod_full,          # logit_algae ~ Warming*Pollution ...
  barn_mod_full_orig,      # barnacles uses prop_algae predictor (mixed scale)
  graz_mod_full_orig
)

print(summary(sem_full))
print(fisherC(sem_full))

claims_full <- dSep(sem_full)
claims_full <- claims_full[order(claims_full$P.Value), ]
print(claims_full)

# 5.2 Reduced + harmonised algae scale everywhere
dat_algae     <- dat_algae     %>% mutate(logit_algae = qlogis(clip01(AvgAlgae / 100)))
dat_barnacles <- dat_barnacles %>% mutate(logit_algae = qlogis(clip01(AvgAlgae / 100)))
dat_grazers   <- dat_grazers   %>% mutate(logit_algae = qlogis(clip01(AvgAlgae / 100)))

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

sem_red <- psem(algae_mod_red, barn_mod_red, graz_mod_logit)

print(summary(sem_red))
print(fisherC(sem_red))

claims_red <- dSep(sem_red)
claims_red <- claims_red[order(claims_red$P.Value), ]
print(claims_red)

## 6. Common-site sensitivity (intersection across submodels)
#   Investigate whether unequal site filtering (Chile removed for algae, Argentina removed for barnacles) drives the results 
#   We refit reduced SEM on the intersection of sites.

common_sites <- Reduce(
  base::intersect,
  list(
    unique(as.character(dat_algae$site)),
    unique(as.character(dat_barnacles$site)),
    unique(as.character(dat_grazers$site))
  )
)

datA_int <- subset(dat_algae,     site %in% common_sites)
datB_int <- subset(dat_barnacles, site %in% common_sites)
datG_int <- subset(dat_grazers,   site %in% common_sites)

algae_mod_int <- lmer(
  logit_algae ~ Warming + Pollution + (1 + Pollution | site),
  data = datA_int
)
barn_mod_int <- lmer(
  log_AvgBarnacles ~ Warming + Pollution + logit_algae + (1 + Pollution | site),
  data = datB_int
)
graz_mod_int <- lmer(
  log_AvgGrazers ~ logit_algae + Pollution + (1 + Pollution | site),
  data = datG_int
)

sem_common <- psem(algae_mod_int, barn_mod_int, graz_mod_int)

print(summary(sem_common))
print(fisherC(sem_common))

## 7. Alternative SEM specifications (additional sensitivities)

# 7.1 RI-only (no random slopes); evidence that it worsens fit / inflates FE
algae_mod_RI <- lmer(
  logit_algae ~ Warming + Pollution + (1 | site),
  data = dat_algae
)

barn_mod_RI <- lmer(
  log_AvgBarnacles ~ Warming + Pollution + logit_algae + (1 | site),
  data = dat_barnacles
)

graz_mod_RI <- lmer(
  log_AvgGrazers ~ logit_algae + Pollution + (1 | site),
  data = dat_grazers
)

sem_RI <- psem(algae_mod_RI, barn_mod_RI, graz_mod_RI)

print(summary(sem_RI))
print(fisherC(sem_RI))
print(AIC(sem_red, sem_RI))

# 7.2 Competition as residual correlation: Algae %~~% Barnacles
# Alternative to directed Algae -> Barnacles
barn_mod_corr <- lmer(
  log_AvgBarnacles ~ Warming + Pollution + (1 + Pollution | site),
  data = dat_barnacles
)

sem_compCorr <- psem(
  algae_mod_red,
  barn_mod_corr,
  graz_mod_logit,
  logit_algae %~~% log_AvgBarnacles
)

print(summary(sem_compCorr))
print(fisherC(sem_compCorr))
print(AIC(sem_red, sem_compCorr))

# 7.3 Barnacles %~~% Grazers residual correlation
# This is retained in the main script
sem_BGcorr <- psem(
  algae_mod_red,
  barn_mod_red,
  graz_mod_logit,
  log_AvgBarnacles %~~% log_AvgGrazers
)

print(summary(sem_BGcorr))
print(fisherC(sem_BGcorr))

claims_bg <- dSep(sem_BGcorr)
claims_bg <- claims_bg[order(claims_bg$P.Value), ]
print(claims_bg)

## 8. Diagrammatic SEM visualisation

# Coefficients from sem_BGcorr
co_raw <- piecewiseSEM::coefs(sem_BGcorr, standardize = "scale")
co <- as_tibble(co_raw, .name_repair = "unique")

pred_col <- names(co)[grepl("^Pred", names(co), ignore.case = TRUE)][1]
resp_col <- names(co)[grepl("^Resp", names(co), ignore.case = TRUE)][1]
p_col    <- names(co)[grepl("^P(\\.|_|$)|P.value", names(co), ignore.case = TRUE)][1]
std_col  <- names(co)[grepl("Std(\\.|_)?Est", names(co), ignore.case = TRUE)][1]
if (is.na(std_col) || std_col == "") {
  std_col <- names(co)[grepl("std", names(co), ignore.case = TRUE) &
                         !grepl("error|se", names(co), ignore.case = TRUE)][1]
}

vars <- c("Warming", "Pollution", "logit_algae", "log_AvgBarnacles", "log_AvgGrazers")
norm <- function(x) gsub("\\d+$", "", x)

paths <- co %>%
  filter(
    !is.na(.data[[pred_col]]),
    .data[[pred_col]] != "(Intercept)",
    norm(.data[[pred_col]]) %in% vars,
    norm(.data[[resp_col]]) %in% vars
  ) %>%
  transmute(
    from = norm(.data[[pred_col]]),
    to   = norm(.data[[resp_col]]),
    std  = suppressWarnings(as.numeric(.data[[std_col]])),
    p    = suppressWarnings(as.numeric(.data[[p_col]]))
  ) %>%
  mutate(
    signif    = p < 0.05,
    color     = ifelse(signif, "black", "gray60"),
    lab_plain = sprintf("%.2f", std),
    lab_html  = ifelse(signif, lab_plain, paste0("<I>", lab_plain, "</I>"))
  )

edges_txt <- glue_collapse(
  glue('  "{paths$from}" -> "{paths$to}" [label=<{paths$lab_html}>, color="{paths$color}"];'),
  sep = "\n"
)

node_labels <- glue_collapse(glue(
  '  "{c("Warming","Pollution","logit_algae","log_AvgBarnacles","log_AvgGrazers")}" [label="{c(
    "Warming",
    "Pollution",
    "Algae coverage (logit)",
    "Mean barnacle count (log)",
    "Mean grazer count (log)"
  )}"];'
), sep = "\n")

cov_txt <- '  "log_AvgBarnacles" -> "log_AvgGrazers" [
      dir=both, style="dashed", color="black",
      label="Residual covariance", constraint=false
  ];'

dot <- glue('
digraph sem {{
  graph [rankdir=TB, splines=true, nodesep=0.6, ranksep=0.8, bgcolor="white", margin=0.1];
  node  [shape=box, style="rounded,filled", fillcolor="white", color="black",
         fontname="Helvetica", fontsize=12, margin="0.06,0.04"];
  edge  [fontname="Helvetica", fontsize=10, arrowsize=0.8, color="black"];

{node_labels}

  {{ rank=same; "Warming"; "Pollution"; }}
  {{ rank=same; "logit_algae"; }}
  {{ rank=same; "log_AvgBarnacles"; "log_AvgGrazers"; }}

{edges_txt}
{cov_txt}
}}')

grViz(dot)