# ============================================================================
# R/15_placebo_period.R   (v2: robust extraction, no pvalue() dependency)
# ----------------------------------------------------------------------------
# Placebo period test: fake treatment dates 2017/2018/2019 + REAL 2020.
#
# FIX vs v1: extracts beta/SE/p from summary(m)$coeftable, avoiding the
# pvalue() function that conflicts with scales::round_any in some R sessions.
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("data.table", "fixest")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(data.table)
library(fixest)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

# Helper: robust extraction
get_estimate <- function(m, var) {
  ct <- summary(m)$coeftable
  if (!var %in% rownames(ct)) {
    return(list(beta = NA, se = NA, p = NA))
  }
  list(
    beta = ct[var, "Estimate"],
    se   = ct[var, "Std. Error"],
    p    = ct[var, "Pr(>|t|)"]
  )
}

# ---- 1. Load panel ---------------------------------------------------------

panel <- readRDS("data/processed/panel_county_year.rds")
panel <- as.data.table(panel)

state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
ia_panel <- panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia")]
stopifnot(nrow(ia_panel) > 0)

cat(sprintf("Iowa panel: %d county-year obs\n", nrow(ia_panel)))

# ---- 2. Placebo loop ------------------------------------------------------

placebo_years <- c(2017, 2018, 2019)
results <- data.table()

cat("\n=== Placebo period test ===\n")
cat("Outcome: log(drug arrests + 1), Iowa-only DiD\n\n")

for (placebo_year in placebo_years) {
  
  # Use only pre-2020 data so the fake date is the only treatment in the window
  sub <- ia_panel[year < 2020]
  
  sub[, fake_post := as.integer(year >= placebo_year)]
  sub[, fake_treat := border * fake_post]
  
  if (sum(sub$fake_post == 1) == 0 || sum(sub$fake_post == 0) == 0) {
    cat(sprintf("  Skipping %d: no variation\n", placebo_year))
    next
  }
  
  m <- feols(
    ln_drug ~ fake_treat + ln_pop + ln_inc | county_fips + year,
    data = sub,
    cluster = ~county_fips
  )
  
  est <- get_estimate(m, "fake_treat")
  
  results <- rbind(results, data.table(
    placebo_year = placebo_year,
    sample = paste0("2015-", placebo_year - 1, " vs ", placebo_year, "-2019"),
    beta = round(est$beta, 4),
    se   = round(est$se, 4),
    p    = round(est$p, 4),
    n    = nobs(m),
    significant_5pct = !is.na(est$p) && est$p < 0.05
  ))
  
  cat(sprintf("  Placebo year %d:  beta = %+.4f  SE = %.4f  p = %.4f  %s\n",
              placebo_year, est$beta, est$se, est$p,
              ifelse(!is.na(est$p) && est$p < 0.05, "*** SIG ***", "(ns)")))
}

# Real 2020 for comparison
if (!"border_x_post" %in% names(ia_panel)) {
  ia_panel[, border_x_post := border * post]
}
m_real <- feols(
  ln_drug ~ border_x_post + ln_pop + ln_inc | county_fips + year,
  data = ia_panel,
  cluster = ~county_fips
)
est_real <- get_estimate(m_real, "border_x_post")

results <- rbind(results, data.table(
  placebo_year = 2020,
  sample = "2015-2022 (REAL, full panel)",
  beta = round(est_real$beta, 4),
  se   = round(est_real$se, 4),
  p    = round(est_real$p, 4),
  n    = nobs(m_real),
  significant_5pct = !is.na(est_real$p) && est_real$p < 0.05
), fill = TRUE)

cat(sprintf("\n  REAL 2020:       beta = %+.4f  SE = %.4f  p = %.4f  %s\n",
            est_real$beta, est_real$se, est_real$p,
            ifelse(!is.na(est_real$p) && est_real$p < 0.05, "*** SIG (expected) ***", "?? not sig ??")))

# ---- 3. Save outputs -------------------------------------------------------

saveRDS(results, "out/placebo_period.rds")

sink("out/placebo_period.txt")
cat("Placebo period test: fake treatment dates in Iowa DiD\n")
cat("=====================================================\n\n")
print(results)
cat("\nInterpretation:\n")
cat(" - If placebo years 2017, 2018, 2019 yield non-significant beta:\n")
cat("   → The 2020 effect is specific to that year, not a pre-existing\n")
cat("     trend or omitted variable.\n")
cat(" - If any placebo year is significant:\n")
cat("   → Pre-existing trend that the DiD captures, weakening the causal\n")
cat("     interpretation.\n")
sink()

cat("\n=== Saved ===\n")
cat("  out/placebo_period.rds\n")
cat("  out/placebo_period.txt\n")

# ---- 4. Defense one-liner --------------------------------------------------

placebo_results <- results[placebo_year < 2020]
n_placebo_sig <- sum(placebo_results$significant_5pct, na.rm = TRUE)

cat("\n", strrep("-", 70), "\n", sep = "")
cat("DEFENSE ONE-LINER:\n")
cat(sprintf("'I run placebo period tests with fake treatment dates 2017,\n"))
cat(sprintf(" 2018, and 2019. %d of 3 placebo estimates are significant at\n", n_placebo_sig))
cat(sprintf(" the 5%% level.\n"))
if (n_placebo_sig == 0) {
  cat(" The 2020 treatment effect is specific to that year, not driven by\n")
  cat(" pre-existing trends. This directly addresses the omitted-variable\n")
  cat(" concern.'\n")
} else {
  cat(" significant placebo estimates indicate pre-existing trends, which\n")
  cat(" the thesis discusses as a limitation.'\n")
}
cat(strrep("-", 70), "\n", sep = "")

cat("\nDone. Add as Section 5.8 (new) in thesis Ch5.\n")
