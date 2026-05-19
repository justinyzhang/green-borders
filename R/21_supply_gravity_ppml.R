# =============================================================================
# 21_supply_gravity_ppml.R
# =============================================================================
# Chapter 4 supply-side gravity regression on Illinois adult-use dispensary
# license entry. Implements the Lychagin directive (5/14 meeting, cluster 6):
#
#   Entry_c = alpha + beta * Border100_c
#           + gamma_IL  * log(PopIL_{c,r}   + 1)
#           + gamma_PROH * log(PopPROH_{c,r} + 1)
#           + delta * X_c + eps_c
#
# Estimated by PPML (Santos Silva & Tenreyro 2006) using fixest::fepois.
# Structural prediction: gamma_PROH > gamma_IL > 0. Cross-border-demand
# coefficient is the differential gamma_PROH - gamma_IL.
#
# Reports three control samples in three columns (Lychagin "sep regs every reg"):
#   (1) All 102 IL counties (baseline)
#   (2) IL counties bordering any state (border-only restriction)
#   (3) Cross-state placebo: counties in IA, MO (pre-2023), IN, KY, WI bordering
#       OTHER prohibition states. Prediction: gamma_PROH null in this sample.
#
# Inputs:
#   - data/il_entry_data.csv (existing: 46 IL border counties with license counts;
#     full 102-county data should be available in green-borders/data/)
#   - out/catchment_pop_il_focal.csv (from 20_catchment_populations.R)
#   - out/catchment_pop_placebo_focal.csv (from 20_catchment_populations.R)
#
# Outputs:
#   - out/gravity_ppml_results.csv (point estimates, SE, p, sample size by spec)
#   - out/gravity_table.tex        (booktabs-formatted LaTeX table for the paper)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(modelsummary)
})

OUT_DIR <- "out"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---- 1. Load entry data ----------------------------------------------------

# Existing IL entry dataset. Adjust path to your repo structure.
# Expected columns: GEOID, COUNAME, licenses (count 2019-2024), pop_2020, etc.
il_entry <- read_csv("data/il_entry_data.csv")

# Catchments from Script 20
cat_il <- read_csv("out/catchment_pop_il_focal.csv")
cat_placebo <- read_csv("out/catchment_pop_placebo_focal.csv")

# ---- 2. Reshape catchment to wide form (one row per county) ----------------

cat_il_wide <- cat_il %>%
  pivot_wider(
    id_cols = c(GEOID, COUNAME, STNAME, county_pop_2020),
    names_from = radius_mi,
    values_from = c(pop_IL, pop_PROH),
    names_glue = "{.value}_{radius_mi}mi"
  )

cat_placebo_wide <- cat_placebo %>%
  pivot_wider(
    id_cols = c(GEOID, COUNAME, STNAME, county_pop_2020),
    names_from = radius_mi,
    values_from = pop_PROH_NEIGHBORS,
    names_glue = "popPROH_{radius_mi}mi"
  )

# ---- 3. Build IL analysis dataset ------------------------------------------

# Merge catchments into entry data; add border indicator at 100 mi
il_data <- il_entry %>%
  left_join(cat_il_wide, by = "GEOID") %>%
  mutate(
    border100   = as.integer(pop_PROH_100mi > 0),
    border50    = as.integer(pop_PROH_50mi > 0),
    border15    = as.integer(pop_PROH_15mi > 0),
    log_popIL_100   = log(pop_IL_100mi   + 1),
    log_popIL_50    = log(pop_IL_50mi    + 1),
    log_popIL_15    = log(pop_IL_15mi    + 1),
    log_popPROH_100 = log(pop_PROH_100mi + 1),
    log_popPROH_50  = log(pop_PROH_50mi  + 1),
    log_popPROH_15  = log(pop_PROH_15mi  + 1),
    log_county_pop  = log(county_pop_2020 + 1)
  )

cat("IL analysis sample: ", nrow(il_data), " counties\n", sep = "")
cat("  with prohibition catchment > 0 at 100 mi: ",
    sum(il_data$border100 == 1), "\n", sep = "")
cat("  with prohibition catchment > 0 at  15 mi: ",
    sum(il_data$border15 == 1), "\n", sep = "")

# ---- 4. Specifications ------------------------------------------------------

# Column 1: all 102 IL counties, three radii
m1_15  <- fepois(licenses ~ log_popIL_15  + log_popPROH_15  + log_county_pop,
                 data = il_data, vcov = "hetero")
m1_50  <- fepois(licenses ~ log_popIL_50  + log_popPROH_50  + log_county_pop,
                 data = il_data, vcov = "hetero")
m1_100 <- fepois(licenses ~ log_popIL_100 + log_popPROH_100 + log_county_pop,
                 data = il_data, vcov = "hetero")

# Column 2: 46 border-only IL counties (border at 100 mi)
il_border_only <- il_data %>% filter(border100 == 1)
m2_15  <- fepois(licenses ~ log_popIL_15  + log_popPROH_15  + log_county_pop,
                 data = il_border_only, vcov = "hetero")
m2_50  <- fepois(licenses ~ log_popIL_50  + log_popPROH_50  + log_county_pop,
                 data = il_border_only, vcov = "hetero")
m2_100 <- fepois(licenses ~ log_popIL_100 + log_popPROH_100 + log_county_pop,
                 data = il_border_only, vcov = "hetero")

# Column 3: PRE-2020 MEDICAL placebo within Illinois. This is the cleanest
# falsification that requires no external data. The IL Medical Cannabis Pilot
# Program issued dispensary licenses 2015-2019; if the gravity coefficient
# gamma_PROH is large and positive on PRE-2020 medical licenses, the
# post-2020 recreational result is contaminated by a pre-existing spatial
# pattern (prohibition borders attract retail activity for reasons unrelated
# to legalization). If gamma_PROH is null on pre-2020 medical, the post-2020
# recreational result is identified by the legalization shock.
#
# Required input: il_entry must have a column `medical_licenses_pre2020`
# (count of MCPP dispensary licenses issued before Jan 1 2020). If absent,
# we attempt to derive it from any `first_license_date` field; failing that,
# we skip the column and print a clear message for Justin to add the data.
#
# Optional secondary placebo: cross-state placebo using counties in non-IL
# prohibition states (IA, MO pre-2023, IN, KY, WI) bordering OTHER
# prohibition states. Requires entry counts in those states (medical
# dispensaries in MO, IA medical CBD locations, etc.). Loaded from
# data/placebo_entry_data.csv if available.

# 3a. Pre-2020 medical placebo (preferred)
medical_col <- NULL
if ("medical_licenses_pre2020" %in% names(il_entry)) {
  medical_col <- "medical_licenses_pre2020"
} else if ("first_license_date" %in% names(il_entry)) {
  # Try to derive from license-level dates if il_entry is dispensary-level
  message("[Attempting to derive medical_licenses_pre2020 from first_license_date]")
  med_counts <- il_entry %>%
    filter(as.Date(first_license_date) < as.Date("2020-01-01")) %>%
    count(GEOID, name = "medical_licenses_pre2020")
  il_data <- il_data %>%
    left_join(med_counts, by = "GEOID") %>%
    mutate(medical_licenses_pre2020 = coalesce(medical_licenses_pre2020, 0L))
  medical_col <- "medical_licenses_pre2020"
}

if (!is.null(medical_col)) {
  cat("\nPlacebo (pre-2020 medical) license counts by border status:\n")
  print(il_data %>%
          group_by(border100) %>%
          summarise(n_counties = n(),
                    total_med_licenses = sum(.data[[medical_col]], na.rm = TRUE),
                    .groups = "drop"))

  m3_15  <- fepois(as.formula(paste(medical_col,
                                    "~ log_popIL_15  + log_popPROH_15  + log_county_pop")),
                   data = il_data, vcov = "hetero")
  m3_50  <- fepois(as.formula(paste(medical_col,
                                    "~ log_popIL_50  + log_popPROH_50  + log_county_pop")),
                   data = il_data, vcov = "hetero")
  m3_100 <- fepois(as.formula(paste(medical_col,
                                    "~ log_popIL_100 + log_popPROH_100 + log_county_pop")),
                   data = il_data, vcov = "hetero")
} else {
  message("[medical_licenses_pre2020 column not found and cannot derive; ",
          "pre-2020 medical placebo column skipped. ",
          "Add this column to data/il_entry_data.csv to enable column 3.]")
  m3_15 <- m3_50 <- m3_100 <- NULL
}

# 3b. Cross-state placebo (secondary, optional)
m4_15 <- m4_50 <- m4_100 <- NULL
placebo_entry_path <- "data/placebo_entry_data.csv"
if (file.exists(placebo_entry_path)) {
  placebo_entry <- read_csv(placebo_entry_path, show_col_types = FALSE)
  placebo_data <- placebo_entry %>%
    left_join(cat_placebo_wide, by = "GEOID") %>%
    mutate(
      log_popPROH_100 = log(popPROH_100mi + 1),
      log_popPROH_50  = log(popPROH_50mi  + 1),
      log_popPROH_15  = log(popPROH_15mi  + 1),
      log_county_pop  = log(county_pop_2020 + 1)
    )
  m4_100 <- fepois(licenses ~ log_popPROH_100 + log_county_pop,
                   data = placebo_data, vcov = "hetero")
  m4_50  <- fepois(licenses ~ log_popPROH_50  + log_county_pop,
                   data = placebo_data, vcov = "hetero")
  m4_15  <- fepois(licenses ~ log_popPROH_15  + log_county_pop,
                   data = placebo_data, vcov = "hetero")
}

# ---- 5. Test gamma_PROH > gamma_IL -----------------------------------------

cat("\n=== Structural test: gamma_PROH > gamma_IL ===\n")
test_differential <- function(model, label) {
  coefs <- coef(model)
  vcv   <- vcov(model)
  i_proh <- grep("popPROH", names(coefs))
  i_il   <- grep("popIL",   names(coefs))
  if (length(i_proh) == 0 || length(i_il) == 0) return(NULL)
  diff <- coefs[i_proh] - coefs[i_il]
  se   <- sqrt(vcv[i_proh, i_proh] + vcv[i_il, i_il] - 2 * vcv[i_proh, i_il])
  t    <- diff / se
  p    <- 2 * pnorm(-abs(t))  # two-sided
  p_1s <- pnorm(-t)            # one-sided: H1: gamma_PROH > gamma_IL
  cat(sprintf("  %s:  diff = %+0.3f  SE = %0.3f  t = %+0.2f  p_1s = %0.3f\n",
              label, diff, se, t, p_1s))
  tibble(spec = label, diff = diff, se = se, t = t, p_1sided = p_1s)
}
diff_tests <- bind_rows(
  test_differential(m1_15,  "All IL, r=15"),
  test_differential(m1_50,  "All IL, r=50"),
  test_differential(m1_100, "All IL, r=100"),
  test_differential(m2_15,  "Border IL, r=15"),
  test_differential(m2_50,  "Border IL, r=50"),
  test_differential(m2_100, "Border IL, r=100")
)

# ---- 6. Compile results ----------------------------------------------------

models <- list(
  "All IL, r=15"     = m1_15,  "All IL, r=50"     = m1_50,  "All IL, r=100"   = m1_100,
  "Border IL, r=15"  = m2_15,  "Border IL, r=50"  = m2_50,  "Border IL, r=100" = m2_100
)
if (!is.null(m3_100)) {
  models <- c(models, list("MedPlacebo, r=15"  = m3_15,
                           "MedPlacebo, r=50"  = m3_50,
                           "MedPlacebo, r=100" = m3_100))
}
if (!is.null(m4_100)) {
  models <- c(models, list("CrossPlacebo, r=15"  = m4_15,
                           "CrossPlacebo, r=50"  = m4_50,
                           "CrossPlacebo, r=100" = m4_100))
}

extract_results <- function(m, label) {
  if (is.null(m)) return(NULL)
  s <- summary(m)
  coefs <- s$coeftable
  tibble(
    spec  = label,
    term  = rownames(coefs),
    beta  = coefs[, "Estimate"],
    se    = coefs[, "Std. Error"],
    z     = coefs[, "z value"],
    p     = coefs[, "Pr(>|z|)"],
    N     = s$nobs
  )
}
results <- imap_dfr(models, ~ extract_results(.x, .y))
write_csv(results, file.path(OUT_DIR, "gravity_ppml_results.csv"))

# ---- 7. LaTeX table --------------------------------------------------------

# modelsummary booktabs export
ms_table <- modelsummary(
  models,
  output = "latex",
  stars = c("*" = 0.10, "**" = 0.05, "***" = 0.01),
  coef_map = c(
    "log_popIL_15"   = "log Pop IL (15 mi)",
    "log_popPROH_15" = "log Pop PROH (15 mi)",
    "log_popIL_50"   = "log Pop IL (50 mi)",
    "log_popPROH_50" = "log Pop PROH (50 mi)",
    "log_popIL_100"   = "log Pop IL (100 mi)",
    "log_popPROH_100" = "log Pop PROH (100 mi)",
    "log_county_pop" = "log County pop"
  ),
  gof_omit = "IC|Log|RMSE|Pseudo R2 Adj.|R2",
  fmt = 3,
  notes = paste(
    "Notes: PPML estimates of dispensary license entry on log-population catchments.",
    "Cross-border-demand structural prediction: gamma_PROH > gamma_IL.",
    "All IL = 102 counties (entry of 2019-2024 recreational licenses).",
    "Border IL = subsample with positive prohibition catchment within 100 mi.",
    "MedPlacebo = same 102 IL counties, outcome is PRE-2020 medical dispensary",
    "licenses (MCPP, 2015-2019). Null prediction: gamma_PROH approximately zero",
    "because medical licenses did not respond to cross-border recreational demand.",
    "CrossPlacebo = counties in IA, MO (pre-2023), IN, KY, WI bordering other",
    "prohibition states (requires data/placebo_entry_data.csv).",
    "Heteroskedasticity-robust standard errors throughout."
  )
)
writeLines(ms_table, file.path(OUT_DIR, "gravity_table.tex"))

write_csv(diff_tests, file.path(OUT_DIR, "gravity_differential_tests.csv"))

# ---- 8. Missouri Feb-2023 within-IL robustness (Tier 2) --------------------

# Restrict to IL counties bordering Missouri. Compare license entry rate
# pre-Feb-2023 vs post-Feb-2023. Prediction: entry growth slows after MO
# legalizes because cross-border demand from MO collapses.
#
# This requires license data with opening-date timestamps. If your existing
# il_entry data has license-level dates, build a county-month panel of entry
# counts and run a within-IL DiD here.

if ("first_license_date" %in% names(il_entry)) {
  message("[Missouri 2023 robustness: would require county-month panel; ",
          "see Phase 2 work plan.]")
}

cat("\n[DONE] gravity_ppml_results.csv, gravity_differential_tests.csv, ",
    "gravity_table.tex written to ", OUT_DIR, "/\n", sep = "")
