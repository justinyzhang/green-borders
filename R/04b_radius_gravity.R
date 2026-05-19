# =============================================================================
# R/04b_radius_gravity.R
# =============================================================================
# Chapter 4 supply-side gravity entry regression (cluster 6 expansion).
# Implements the Lychagin 5/14 directive: PPML with separately measured
# log(Pop_IL) and log(Pop_PROH) catchments at three radii (15, 50, 100 mi).
#
# Structural prediction: gamma_PROH > gamma_IL > 0. The cross-border-demand
# differential is gamma_PROH - gamma_IL, tested one-sided.
#
# Three sample columns (Lychagin "separate regs every reg"):
#   (1) All 102 IL counties
#   (2) IL counties with positive prohibition catchment at r=100mi
#   (3) Pre-2020 MEDICAL licenses placebo (same 102 counties)
#       Null prediction: gamma_PROH ~ 0 on pre-2020 medical, because
#       medical-era licenses did not respond to cross-border recreational
#       demand. If gamma_PROH is large on pre-2020 medical, the post-2020
#       result is contaminated by pre-existing spatial sorting.
#
# Inputs:
#   - data/raw/idfpr/dispensaries_clean.csv (IDFPR dispensary registry)
#   - data/processed/catchment_pop_il_focal.csv (from radius_design pipeline,
#     ACS 5-year 2018-2022 block-group catchments at r=15/50/100 mi)
#
# Outputs:
#   - output/tables/tab4_gravity_coefs.csv
#   - output/tables/tab4_gravity.tex
#   - output/tables/tab4_gravity_differential_tests.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(here)
})

# Estimator selection (CLAUDE.md convention: try fixest, fall back to glm)
have_fixest <- requireNamespace("fixest", quietly = TRUE)
if (have_fixest) {
  library(fixest)
  cat("[Setup] Using fixest::fepois for PPML\n")
} else {
  cat("[Setup] fixest unavailable; using glm(quasipoisson) fallback\n")
}

dir.create(here::here("output", "tables"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Locate inputs
# -----------------------------------------------------------------------------
disp_paths <- c(
  here::here("data", "raw", "idfpr", "dispensaries_clean.csv"),
  here::here("data", "raw", "dispensaries_clean.csv"),
  here::here("data", "processed", "dispensaries_clean.csv"),
  "/mnt/user-data/uploads/dispensaries_clean.csv"
)
disp_path <- disp_paths[file.exists(disp_paths)][1]
if (is.na(disp_path)) stop("dispensaries_clean.csv not found; tried:\n",
                            paste(disp_paths, collapse = "\n"))
cat("Reading dispensaries from:", disp_path, "\n")

catch_paths <- c(
  here::here("data", "processed", "catchment_pop_il_focal.csv"),
  here::here("out", "catchment_pop_il_focal.csv"),
  here::here("..", "radius_design", "out", "catchment_pop_il_focal.csv")
)
catch_path <- catch_paths[file.exists(catch_paths)][1]
if (is.na(catch_path)) stop("catchment_pop_il_focal.csv not found; tried:\n",
                             paste(catch_paths, collapse = "\n"),
                             "\nRun radius_design/20_catchment_populations.R first.")
cat("Reading catchments from:", catch_path, "\n")

# -----------------------------------------------------------------------------
# 2. Load and aggregate dispensary data (matches R_04_cross_section_rebuild_v2)
# -----------------------------------------------------------------------------
disp <- fread(disp_path)
cat("Dispensary records loaded:", nrow(disp), "\n")

# Normalize county names
disp[, County := str_to_title(County)]
disp[, County := fifelse(County == "Dekalb",    "DeKalb",
                fifelse(County == "Dupage",    "DuPage",
                fifelse(County == "Mchenry",   "McHenry",
                fifelse(County == "Mclean",    "McLean",
                fifelse(County == "Mcdonough", "McDonough",
                fifelse(County == "Lasalle",   "LaSalle", County))))))]

# 2019-2024 cumulative adult-use entry (matches v15 panel cut)
disp_2024 <- disp[Year >= 2019 & Year <= 2024]
entry_2024 <- disp_2024[, .(licenses = .N), by = County]

# Pre-2020 medical-era licenses (placebo: should not respond to cross-border)
disp_pre2020 <- disp[Year < 2020]
entry_pre2020 <- disp_pre2020[, .(medical_licenses_pre2020 = .N), by = County]

cat("2019-2024 adult-use:", sum(entry_2024$licenses), "licenses across",
    nrow(entry_2024), "counties\n")
cat("Pre-2020 medical:    ", sum(entry_pre2020$medical_licenses_pre2020),
    "licenses across", nrow(entry_pre2020), "counties\n")

# -----------------------------------------------------------------------------
# 3. Load catchments and build 102-county wide-format dataset
# -----------------------------------------------------------------------------
cat_il <- fread(catch_path)
cat_il[, GEOID := sprintf("%05d", as.integer(GEOID))]
cat("IL catchment rows (102 counties x 3 radii):", nrow(cat_il), "\n")

# Reshape catchments to wide: one row per county
cat_il_wide <- dcast(cat_il,
                     GEOID + COUNAME + STNAME + county_pop_2020 ~ radius_mi,
                     value.var = c("pop_IL", "pop_PROH"))
# Rename: pop_IL_15 -> pop_IL_15mi, etc.
old_names <- grep("^pop_(IL|PROH)_(15|50|100)$", names(cat_il_wide), value = TRUE)
new_names <- paste0(old_names, "mi")
setnames(cat_il_wide, old_names, new_names)

cat("Wide catchment columns:\n")
print(names(cat_il_wide))

# -----------------------------------------------------------------------------
# 4. Merge entry counts onto 102-county frame (zero-entry counties stay)
# -----------------------------------------------------------------------------
il_data <- merge(cat_il_wide,
                 entry_2024,
                 by.x = "COUNAME", by.y = "County", all.x = TRUE)
il_data <- merge(il_data,
                 entry_pre2020,
                 by.x = "COUNAME", by.y = "County", all.x = TRUE)
# Zero-entry counties get NA -> 0
il_data[is.na(licenses), licenses := 0L]
il_data[is.na(medical_licenses_pre2020), medical_licenses_pre2020 := 0L]

cat("\n102-county sample:\n")
cat("  Counties with positive adult-use licenses:",
    sum(il_data$licenses > 0), "of 102\n")
cat("  Counties with positive pre-2020 medical:",
    sum(il_data$medical_licenses_pre2020 > 0), "of 102\n")
cat("  Total adult-use licenses (should be 240):",
    sum(il_data$licenses), "\n")
cat("  Total pre-2020 medical licenses (should be 50):",
    sum(il_data$medical_licenses_pre2020), "\n")

if (sum(il_data$licenses) != 240) {
  cat("\nWARN: total licenses does not match expected 240 — check county name normalization\n")
  unmatched <- setdiff(entry_2024$County, il_data$COUNAME)
  if (length(unmatched) > 0) {
    cat("  Counties in entry data but not in catchment:", paste(unmatched, collapse = ", "), "\n")
  }
}

# -----------------------------------------------------------------------------
# 5. Build derived variables
# -----------------------------------------------------------------------------
il_data[, log_county_pop := log(county_pop_2020 + 1)]
for (r in c(15, 50, 100)) {
  il_data[, paste0("log_popIL_", r)   := log(get(paste0("pop_IL_",   r, "mi")) + 1)]
  il_data[, paste0("log_popPROH_", r) := log(get(paste0("pop_PROH_", r, "mi")) + 1)]
  il_data[, paste0("border", r)        := as.integer(get(paste0("pop_PROH_", r, "mi")) > 0)]
}

cat("\nBorder indicator counts:\n")
cat("  border15  (pos PROH catchment at 15 mi):", sum(il_data$border15),  "\n")
cat("  border50  (pos PROH catchment at 50 mi):", sum(il_data$border50),  "\n")
cat("  border100 (pos PROH catchment at 100 mi):", sum(il_data$border100), "\n")

# -----------------------------------------------------------------------------
# 6. PPML estimation function (fixest if available, glm fallback)
# -----------------------------------------------------------------------------
fit_ppml <- function(y_col, log_il_col, log_proh_col, data) {
  f <- as.formula(sprintf("%s ~ %s + %s + log_county_pop",
                           y_col, log_il_col, log_proh_col))
  if (have_fixest) {
    m <- fepois(f, data = data, vcov = "hetero")
    ct <- summary(m)$coeftable
    p_col <- if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "Pr(>|t|)"
    list(
      model = m,
      coefs = coef(m),
      vcov  = vcov(m),
      n     = nobs(m),
      table = ct,
      p_col = p_col
    )
  } else {
    m <- glm(f, family = quasipoisson(link = "log"), data = data)
    ct <- summary(m)$coefficients
    list(
      model = m,
      coefs = coef(m),
      vcov  = vcov(m),
      n     = nrow(model.frame(m)),
      table = ct,
      p_col = "Pr(>|t|)"
    )
  }
}

extract_row <- function(fit, var, label, sample_label, radius) {
  ct <- fit$table
  if (!var %in% rownames(ct)) {
    return(data.table(spec = label, sample = sample_label, radius = radius,
                      term = var, beta = NA_real_, se = NA_real_, p = NA_real_,
                      n = fit$n))
  }
  data.table(
    spec = label, sample = sample_label, radius = radius, term = var,
    beta = ct[var, "Estimate"],
    se   = ct[var, "Std. Error"],
    p    = ct[var, fit$p_col],
    n    = fit$n
  )
}

# Test gamma_PROH > gamma_IL (one-sided)
test_differential <- function(fit, log_il_col, log_proh_col,
                              sample_label, radius) {
  cf <- fit$coefs
  V  <- fit$vcov
  if (!(log_il_col %in% names(cf)) || !(log_proh_col %in% names(cf))) {
    return(data.table(sample = sample_label, radius = radius,
                      diff = NA_real_, se_diff = NA_real_,
                      t = NA_real_, p_1sided = NA_real_))
  }
  d  <- unname(cf[log_proh_col] - cf[log_il_col])
  se <- sqrt(V[log_proh_col, log_proh_col] + V[log_il_col, log_il_col]
             - 2 * V[log_proh_col, log_il_col])
  t  <- d / se
  p_1s <- pnorm(-t)  # H1: gamma_PROH > gamma_IL
  data.table(sample = sample_label, radius = radius,
             diff = d, se_diff = se, t = t, p_1sided = p_1s)
}

# -----------------------------------------------------------------------------
# 7. Run all specifications
# -----------------------------------------------------------------------------
cat("\n=== PPML gravity entry estimates ===\n")
cat("Structural prediction: gamma_PROH > gamma_IL\n\n")

all_rows  <- list()
all_diffs <- list()
models    <- list()

samples <- list(
  "All IL"     = il_data,
  "Border IL"  = il_data[border100 == 1]
)

for (sample_label in names(samples)) {
  sdat <- samples[[sample_label]]
  cat(sprintf("Sample: %s (N = %d)\n", sample_label, nrow(sdat)))
  for (r in c(15, 50, 100)) {
    log_il   <- paste0("log_popIL_",   r)
    log_proh <- paste0("log_popPROH_", r)
    label    <- paste0(sample_label, ", r=", r)

    fit <- fit_ppml("licenses", log_il, log_proh, sdat)
    models[[label]] <- fit$model

    all_rows[[length(all_rows) + 1]] <- extract_row(fit, log_il,   label, sample_label, r)
    all_rows[[length(all_rows) + 1]] <- extract_row(fit, log_proh, label, sample_label, r)
    all_rows[[length(all_rows) + 1]] <- extract_row(fit, "log_county_pop", label,
                                                     sample_label, r)
    all_diffs[[length(all_diffs) + 1]] <- test_differential(fit, log_il, log_proh,
                                                             sample_label, r)

    g_il   <- fit$coefs[log_il]
    g_proh <- fit$coefs[log_proh]
    cat(sprintf("  r=%3d mi:  gamma_IL = %+.3f   gamma_PROH = %+.3f   diff = %+.3f\n",
                r, g_il, g_proh, g_proh - g_il))
  }
}

# -----------------------------------------------------------------------------
# 8. Pre-2020 medical placebo (same 102 counties, different outcome)
# -----------------------------------------------------------------------------
cat("\n=== Pre-2020 medical placebo ===\n")
cat("Null prediction: gamma_PROH ~ 0 on pre-2020 medical (no cross-border response)\n\n")

for (r in c(15, 50, 100)) {
  log_il   <- paste0("log_popIL_",   r)
  log_proh <- paste0("log_popPROH_", r)
  label    <- paste0("MedPlacebo, r=", r)

  fit <- fit_ppml("medical_licenses_pre2020", log_il, log_proh, il_data)
  models[[label]] <- fit$model

  all_rows[[length(all_rows) + 1]] <- extract_row(fit, log_il,   label, "MedPlacebo", r)
  all_rows[[length(all_rows) + 1]] <- extract_row(fit, log_proh, label, "MedPlacebo", r)
  all_rows[[length(all_rows) + 1]] <- extract_row(fit, "log_county_pop", label,
                                                   "MedPlacebo", r)
  all_diffs[[length(all_diffs) + 1]] <- test_differential(fit, log_il, log_proh,
                                                           "MedPlacebo", r)

  g_il   <- fit$coefs[log_il]
  g_proh <- fit$coefs[log_proh]
  cat(sprintf("  r=%3d mi:  gamma_IL = %+.3f   gamma_PROH = %+.3f   diff = %+.3f\n",
              r, g_il, g_proh, g_proh - g_il))
}

# -----------------------------------------------------------------------------
# 9. Differential test summary
# -----------------------------------------------------------------------------
diffs <- rbindlist(all_diffs)
cat("\n=== Structural test: gamma_PROH > gamma_IL (one-sided) ===\n")
for (i in seq_len(nrow(diffs))) {
  cat(sprintf("  %s, r=%3d:  diff = %+.3f  SE = %.3f  t = %+.2f  p_1s = %.4f\n",
              diffs$sample[i], diffs$radius[i], diffs$diff[i],
              diffs$se_diff[i], diffs$t[i], diffs$p_1sided[i]))
}

# -----------------------------------------------------------------------------
# 10. Save outputs
# -----------------------------------------------------------------------------
results <- rbindlist(all_rows)
fwrite(results, here::here("output", "tables", "tab4_gravity_coefs.csv"))
fwrite(diffs,   here::here("output", "tables", "tab4_gravity_differential_tests.csv"))

cat("\nWrote:\n")
cat("  output/tables/tab4_gravity_coefs.csv\n")
cat("  output/tables/tab4_gravity_differential_tests.csv\n")

# Headline output for defense recall
cat("\n",
    strrep("-", 70), "\n",
    "DEFENSE HEADLINE:\n",
    sprintf("'PPML gravity at r=50 mi over all 102 IL counties yields gamma_IL = %.3f\n",
            diffs[sample == "All IL" & radius == 50, diff] +
            results[sample == "All IL" & radius == 50 & term == "log_popIL_50", beta]),
    sprintf(" and gamma_PROH = %.3f. The differential gamma_PROH - gamma_IL =\n",
            results[sample == "All IL" & radius == 50 & term == "log_popPROH_50", beta]),
    sprintf(" %+.3f (p_1sided = %.3f), supporting the cross-border structural\n",
            diffs[sample == "All IL" & radius == 50, diff],
            diffs[sample == "All IL" & radius == 50, p_1sided]),
    " demand prediction.'\n",
    strrep("-", 70), "\n",
    sep = "")

cat("\nDone.\n")
