# =============================================================================
# R/12b_radius_dose_did.R
# =============================================================================
# Chapter 5 continuous-dose DiD (cluster 11 implementation).
# Replaces the v15 binary IL-border indicator with a CONTINUOUS treatment
# intensity = log(Pop_IL within radius r of county centroid).
#
# Specifications:
#   (A) Continuous-dose DiD at three radii (15, 50, 100 mi)
#   (B) Sanity check: replicate v15 binary spec, beta should be ~0.422
#   (C) Continuous-dose event study with k = -5..+2 at r=50 mi
#
# Implements Lychagin 5/14 directive (cluster 11). Sample: 99 Iowa counties,
# 2015-2022, NIBRS drug arrests + 3 other outcomes.
#
# Uses base lm + sandwich::vcovCL (matches v15 infrastructure; no fixest dep).
#
# Inputs:
#   - data/processed/panel_county_year.rds (from R/03_panel.R)
#   - data/processed/catchment_pop_ia_focal.csv (from radius_design pipeline)
#
# Outputs:
#   - output/tables/tab5_continuous_dose.csv
#   - output/tables/tab5_continuous_dose.tex
#   - output/figures/fig5_continuous_eventstudy.pdf
#   - output/tables/v15_binary_replication.txt  (sanity-check report)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(here)
})

dir.create(here::here("output", "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("output", "figures"), showWarnings = FALSE, recursive = TRUE)

# AER-style plot palette (matches fig3_eventstudy in v15 appendix)
AER_RED  <- "#9E2A2B"
AER_GREY <- "#666666"

# -----------------------------------------------------------------------------
# 1. Load panel and catchments
# -----------------------------------------------------------------------------
panel_path <- here::here("data", "processed", "panel_county_year.rds")
if (!file.exists(panel_path)) {
  panel_path <- "/mnt/user-data/uploads/panel_county_year.rds"
}
stopifnot(file.exists(panel_path))
cat("Reading panel from:", panel_path, "\n")
panel <- as.data.table(readRDS(panel_path))

# Iowa-only sample (cluster 5 spec; state-by-year FE collapses to year FE here)
ia <- panel[state_name == "iowa"]
cat(sprintf("Iowa sample: %d rows, %d counties, %d years\n",
            nrow(ia), uniqueN(ia$county_fips), uniqueN(ia$year)))

catch_paths <- c(
  here::here("data", "processed", "catchment_pop_ia_focal.csv"),
  here::here("out", "catchment_pop_ia_focal.csv"),
  here::here("..", "radius_design", "out", "catchment_pop_ia_focal.csv")
)
catch_path <- catch_paths[file.exists(catch_paths)][1]
if (is.na(catch_path)) stop("catchment_pop_ia_focal.csv not found; tried:\n",
                             paste(catch_paths, collapse = "\n"),
                             "\nRun radius_design/20_catchment_populations.R first.")
cat("Reading catchments from:", catch_path, "\n")

cat_ia <- fread(catch_path)
cat_ia[, GEOID := sprintf("%05d", as.integer(GEOID))]

# Reshape to wide: one row per county with pop_IL at each radius
cat_ia_wide <- dcast(cat_ia,
                     GEOID + COUNAME ~ radius_mi,
                     value.var = c("pop_IL", "pop_IA"))
old_names <- grep("^pop_(IL|IA)_(15|50|100)$", names(cat_ia_wide), value = TRUE)
new_names <- paste0(old_names, "mi")
setnames(cat_ia_wide, old_names, new_names)
setnames(cat_ia_wide, "GEOID", "county_fips")

# Merge into Iowa panel
ia <- merge(ia, cat_ia_wide[, .(county_fips, pop_IL_15mi, pop_IL_50mi, pop_IL_100mi)],
            by = "county_fips", all.x = TRUE)
unmatched <- sum(is.na(ia$pop_IL_50mi))
cat(sprintf("Catchment merge: %d of %d rows matched\n",
            nrow(ia) - unmatched, nrow(ia)))
stopifnot(unmatched == 0)

# -----------------------------------------------------------------------------
# 2. Derive continuous-dose variables
# -----------------------------------------------------------------------------
ia[, dose_15  := log(pop_IL_15mi  + 1)]
ia[, dose_50  := log(pop_IL_50mi  + 1)]
ia[, dose_100 := log(pop_IL_100mi + 1)]
# Center at within-Iowa mean so the dose-by-post coefficient is interpretable
# as the differential effect of one log-unit higher IL exposure.
ia[, dose_15_c  := dose_15  - mean(dose_15,  na.rm = TRUE)]
ia[, dose_50_c  := dose_50  - mean(dose_50,  na.rm = TRUE)]
ia[, dose_100_c := dose_100 - mean(dose_100, na.rm = TRUE)]

# Verify bridge counties have highest doses
bridge_geoids <- c("19045","19057","19061","19097","19111","19115","19139","19163")
cat("\n=== Dose values for v15 bridge counties at r=50 mi ===\n")
print(ia[county_fips %in% bridge_geoids & year == 2019,
         .(county_fips, name, pop_IL_50mi, dose_50)][order(-pop_IL_50mi)])

# -----------------------------------------------------------------------------
# 3. Spec A: Continuous-dose DiD at three radii, four outcomes
# -----------------------------------------------------------------------------
# Spec: ln(Y_ct + 1) = alpha_c + lambda_t + delta * dose_c * Post_t
#                    + b1 * ln_pop_ct + b2 * ln_inc_ct + u_ct
# Cluster SE at county (99 clusters, sandwich HC1).

fit_dose_did <- function(y_col, dose_col, data) {
  data[, dose_x_post := get(dose_col) * post]
  f <- as.formula(sprintf(
    "%s ~ dose_x_post + ln_pop + ln_inc + factor(county_fips) + factor(year)",
    y_col
  ))
  m <- lm(f, data = data)
  vc <- sandwich::vcovCL(m, cluster = data$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  list(
    model = m,
    beta  = unname(ct["dose_x_post", "Estimate"]),
    se    = unname(ct["dose_x_post", "Std. Error"]),
    p     = unname(ct["dose_x_post", "Pr(>|t|)"]),
    ci_lo = unname(ci["dose_x_post", 1]),
    ci_hi = unname(ci["dose_x_post", 2]),
    n     = nobs(m)
  )
}

cat("\n=== Spec A: continuous-dose DiD across 4 outcomes x 3 radii ===\n")
cat("Coefficient is on dose_c x Post_t (interpretation: differential post-2020\n")
cat("change per one log-unit higher IL-side catchment population).\n\n")

outcomes <- c("ln_drug" = "Drug", "ln_owi" = "OWI",
              "ln_property" = "Property", "ln_violent" = "Violent")
radii_doses <- c("dose_15_c" = 15, "dose_50_c" = 50, "dose_100_c" = 100)

dose_rows <- list()
for (out_var in names(outcomes)) {
  for (dose_var in names(radii_doses)) {
    fit <- fit_dose_did(out_var, dose_var, copy(ia))
    dose_rows[[length(dose_rows) + 1]] <- data.table(
      outcome = outcomes[out_var],
      radius_mi = radii_doses[dose_var],
      beta = round(fit$beta, 4),
      se   = round(fit$se, 4),
      p    = round(fit$p, 4),
      ci_lo = round(fit$ci_lo, 4),
      ci_hi = round(fit$ci_hi, 4),
      n    = fit$n
    )
  }
}
dose_results <- rbindlist(dose_rows)
print(dose_results)

fwrite(dose_results,
       here::here("output", "tables", "tab5_continuous_dose.csv"))

# -----------------------------------------------------------------------------
# 4. Spec B: Sanity check - replicate v15 binary specification
# -----------------------------------------------------------------------------
# Should recover beta ~ 0.422, SE ~ 0.198 from v15 Table 5.1 row 5 (drug).

cat("\n=== Spec B: v15 binary replication (sanity check) ===\n")
ia[, treat_v15 := border * post]
m_v15 <- lm(ln_drug ~ treat_v15 + ln_pop + ln_inc
              + factor(county_fips) + factor(year), data = ia)
vc_v15 <- sandwich::vcovCL(m_v15, cluster = ia$county_fips, type = "HC1")
ct_v15 <- lmtest::coeftest(m_v15, vcov. = vc_v15)
beta_v15 <- ct_v15["treat_v15", "Estimate"]
se_v15   <- ct_v15["treat_v15", "Std. Error"]
p_v15    <- ct_v15["treat_v15", "Pr(>|t|)"]

cat(sprintf("v15 binary headline: beta = %.4f, SE = %.4f, p = %.4f\n",
            beta_v15, se_v15, p_v15))
cat(sprintf("v15 thesis Table 5.1: beta = 0.422,  SE = 0.198,  p = 0.033\n"))

agreement <- abs(beta_v15 - 0.422) < 0.05
if (agreement) {
  cat("✓ Replication MATCHES v15 thesis (within 0.05 tolerance)\n")
} else {
  cat(sprintf("⚠ Replication differs by %.3f — investigate panel build\n",
              abs(beta_v15 - 0.422)))
}

sink(here::here("output", "tables", "v15_binary_replication.txt"))
cat("v15 binary headline replication on real panel_county_year.rds\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Replicated:   beta = %.4f, SE = %.4f, p = %.4f, N = %d\n",
            beta_v15, se_v15, p_v15, nobs(m_v15)))
cat(sprintf("v15 reported: beta = 0.422,  SE = 0.198, p = 0.033, N = 792\n"))
cat(sprintf("Match (|diff| < 0.05): %s\n", if (agreement) "YES" else "NO"))
sink()

# Compute implied effect at median bridge dose for headline numeric
bridge_dose_med <- median(ia[county_fips %in% bridge_geoids & year == 2019, dose_50])
interior_dose_med <- median(ia[!(county_fips %in% bridge_geoids) & year == 2019, dose_50])
dose_delta <- bridge_dose_med - interior_dose_med

beta_dose_50_drug <- dose_results[outcome == "Drug" & radius_mi == 50, beta]
implied_v15_equiv <- beta_dose_50_drug * dose_delta

cat(sprintf("\nImplied bridge-vs-interior differential at r=50:\n"))
cat(sprintf("  dose differential (bridge med - interior med) = %.3f\n", dose_delta))
cat(sprintf("  beta_continuous_drug_50 = %.4f\n", beta_dose_50_drug))
cat(sprintf("  implied effect = %.4f x %.3f = %.4f\n",
            beta_dose_50_drug, dose_delta, implied_v15_equiv))
cat(sprintf("  (v15 binary: 0.422; ratio to v15 = %.2f)\n",
            implied_v15_equiv / 0.422))

# -----------------------------------------------------------------------------
# 5. Spec C: Continuous-dose event study, drug at r=50 mi
# -----------------------------------------------------------------------------
cat("\n=== Spec C: continuous-dose event study (drug, r=50 mi) ===\n")

ia[, year_rel := year - 2020]
event_years <- sort(unique(ia$year_rel))
ref_k <- -1

# Build dose-by-event-time interactions (drop k = -1 as reference)
dose_es <- copy(ia)
for (k in event_years) {
  if (k == ref_k) next
  cname <- sprintf("dose50_k%s", if (k < 0) paste0("m", abs(k)) else as.character(k))
  dose_es[, (cname) := as.integer(year_rel == k) * dose_50_c]
}
k_cols <- grep("^dose50_k", names(dose_es), value = TRUE)

f_es <- as.formula(sprintf(
  "ln_drug ~ %s + ln_pop + ln_inc + factor(county_fips) + factor(year)",
  paste(k_cols, collapse = " + ")
))
m_es <- lm(f_es, data = dose_es)
vc_es <- sandwich::vcovCL(m_es, cluster = dose_es$county_fips, type = "HC1")
ct_es <- lmtest::coeftest(m_es, vcov. = vc_es)
ci_es <- lmtest::coefci(m_es, vcov. = vc_es, level = 0.95)

# Tidy
es_tidy <- data.table(
  k = sapply(k_cols, function(s) {
    x <- sub("^dose50_k", "", s)
    if (substr(x, 1, 1) == "m") -as.integer(sub("^m", "", x)) else as.integer(x)
  }),
  est   = ct_es[k_cols, "Estimate"],
  se    = ct_es[k_cols, "Std. Error"],
  ci_lo = ci_es[k_cols, 1],
  ci_hi = ci_es[k_cols, 2]
)
# Add reference row k = -1
es_tidy <- rbindlist(list(es_tidy,
                          data.table(k = -1, est = 0, se = 0,
                                     ci_lo = 0, ci_hi = 0)))
setorder(es_tidy, k)
print(es_tidy)

# Joint Wald test of pre-period coefficients (k = -5..-2)
pre_cols <- grep("^dose50_km[2345]$", names(coef(m_es)), value = TRUE)
wald_pre <- tryCatch({
  R <- matrix(0, nrow = length(pre_cols), ncol = length(coef(m_es)))
  for (j in seq_along(pre_cols)) {
    R[j, which(names(coef(m_es)) == pre_cols[j])] <- 1
  }
  beta_pre <- R %*% coef(m_es)
  V_pre    <- R %*% vc_es %*% t(R)
  W <- as.numeric(t(beta_pre) %*% solve(V_pre) %*% beta_pre)
  F <- W / length(pre_cols)
  p <- pf(F, df1 = length(pre_cols), df2 = nobs(m_es) - length(coef(m_es)),
          lower.tail = FALSE)
  list(F = F, p = p, df1 = length(pre_cols))
}, error = function(e) {
  cat("[Wald test failed:", conditionMessage(e), "]\n")
  NULL
})

if (!is.null(wald_pre)) {
  cat(sprintf("\nJoint Wald test of pre-period dose coefficients (k = -5..-2):\n"))
  cat(sprintf("  F = %.3f, p_joint = %.4f\n", wald_pre$F, wald_pre$p))
}

# -----------------------------------------------------------------------------
# 6. AER-style event study plot
# -----------------------------------------------------------------------------
p_es <- ggplot(es_tidy, aes(x = k, y = est)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted",
             linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, color = AER_GREY, linewidth = 0.5) +
  geom_line(color = AER_RED, linewidth = 0.6) +
  geom_point(color = AER_RED, size = 2.4) +
  scale_x_continuous(breaks = -5:2) +
  labs(
    x = "Years from January 2020 IL legalization",
    y = "Dose-scaled coefficient on log(IL pop within 50 mi)",
    title = NULL
  ) +
  annotate("text", x = -5, y = max(es_tidy$ci_hi) * 0.95,
           label = "Continuous-dose drug arrests, Iowa (r = 50 mi)",
           hjust = 0, size = 4) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

ggsave(here::here("output", "figures", "fig5_continuous_eventstudy.pdf"),
       p_es, width = 7, height = 4.2)
cat("\nWrote output/figures/fig5_continuous_eventstudy.pdf\n")

# Also save event study coefs as CSV
fwrite(es_tidy,
       here::here("output", "tables", "tab5_dose_eventstudy_coefs.csv"))

# -----------------------------------------------------------------------------
# 7. Summary outputs
# -----------------------------------------------------------------------------
cat("\nWrote:\n")
cat("  output/tables/tab5_continuous_dose.csv\n")
cat("  output/tables/v15_binary_replication.txt\n")
cat("  output/tables/tab5_dose_eventstudy_coefs.csv\n")
cat("  output/figures/fig5_continuous_eventstudy.pdf\n")

# Defense one-liner
drug50 <- dose_results[outcome == "Drug" & radius_mi == 50]
cat("\n",
    strrep("-", 70), "\n",
    "DEFENSE ONE-LINER:\n",
    sprintf("'I re-estimate the Iowa DiD with a CONTINUOUS treatment intensity\n"),
    sprintf(" equal to log(Pop_IL within 50 mi). The dose-by-post coefficient on\n"),
    sprintf(" log drug arrests is beta = %.4f (SE %.4f, p = %.3f), with\n",
            drug50$beta, drug50$se, drug50$p),
    sprintf(" identification across all 99 counties rather than 8-vs-91. The\n"),
    sprintf(" v15 binary specification replicates at beta = %.3f (SE %.3f).\n",
            beta_v15, se_v15),
    sprintf(" At median bridge-vs-interior dose, the continuous spec implies\n"),
    sprintf(" %.3f, recovering the v15 binary 0.422 within %.0f%% relative.'\n",
            implied_v15_equiv, 100 * abs(implied_v15_equiv - 0.422) / 0.422),
    strrep("-", 70), "\n",
    sep = "")

cat("\nDone.\n")
