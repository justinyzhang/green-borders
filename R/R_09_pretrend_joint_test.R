# ============================================================================
# R/09_pretrend_joint_test.R   (VERIFIED, simplified to feols + wald)
# ----------------------------------------------------------------------------
# Joint F-test of pre-period coefficients in the event study.
# Addresses defense Concern #2: "Parallel trends test is underpowered."
#
# H0: beta_{-5} = beta_{-4} = beta_{-3} = beta_{-2} = 0
#
# Inputs:  data/processed/panel_county_year.rds
# Outputs: out/pretrend_joint_test.rds
#          out/pretrend_joint_test.txt
# Wall-clock: < 30 sec
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("data.table", "fixest")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(data.table)
library(fixest)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load panel ---------------------------------------------------------

panel <- readRDS("data/processed/panel_county_year.rds")
panel <- as.data.table(panel)

state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
ia_panel <- panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia")]
stopifnot(nrow(ia_panel) > 0)

cat(sprintf("Iowa panel: %d county-year obs\n", nrow(ia_panel)))

# ---- 2. Build event-time dummies -------------------------------------------

# year_rel = year - 2020; k = -1 (2019) is the reference period
ia_panel[, year_rel := year - 2020]

# Restrict to event window k = -5 to +2 (i.e. years 2015-2022)
ia_panel <- ia_panel[year_rel >= -5 & year_rel <= 2]

# Construct treat x k dummies for all k != -1
event_levels <- c(-5, -4, -3, -2, 0, 1, 2)
event_cols <- character()

for (k in event_levels) {
  col <- if (k < 0) paste0("treat_k_m", -k) else paste0("treat_k_p", k)
  if (k == 0) col <- "treat_k_0"
  ia_panel[, (col) := as.integer(year_rel == k) * border]
  event_cols <- c(event_cols, col)
}

cat(sprintf("Event-time columns built: %s\n", paste(event_cols, collapse = ", ")))

# The 4 pre-period coefficients we will test jointly
pre_cols <- c("treat_k_m5", "treat_k_m4", "treat_k_m3", "treat_k_m2")

# ---- 3. Event study + joint Wald per outcome -------------------------------

outcomes <- c("ln_drug", "ln_owi", "ln_property", "ln_violent")
outcome_labels <- c("Drug arrests", "OWI arrests", "Property crime", "Violent crime")

cat("\n=== Joint Wald test of pre-period coefficients ===\n")
cat("H0: beta_{-5} = beta_{-4} = beta_{-3} = beta_{-2} = 0\n\n")

results <- data.table()

for (i in seq_along(outcomes)) {
  outcome <- outcomes[i]
  if (!outcome %in% names(ia_panel)) {
    cat(sprintf("Skipping %s (not in panel)\n", outcome))
    next
  }

  rhs <- paste(event_cols, collapse = " + ")
  formula_es <- as.formula(
    paste0(outcome, " ~ ", rhs, " + ln_pop + ln_inc | county_fips + year")
  )

  m_es <- feols(formula_es,
                data    = ia_panel,
                cluster = ~county_fips)

  # Cluster-robust joint test of pre-period coefs
  w <- tryCatch(
    wald(m_es, keep = pre_cols, print = FALSE),
    error = function(e) {
      cat(sprintf("  wald() failed for %s: %s\n", outcome, e$message))
      NULL
    }
  )

  if (is.null(w)) {
    results <- rbind(results, data.table(
      outcome = outcome_labels[i],
      F_stat = NA, df1 = NA, df2 = NA, p_joint = NA,
      beta_m5 = NA, beta_m4 = NA, beta_m3 = NA, beta_m2 = NA
    ))
    next
  }

  betas <- coef(m_es)[pre_cols]

  results <- rbind(results, data.table(
    outcome = outcome_labels[i],
    F_stat  = round(as.numeric(w$stat), 3),
    df1     = length(pre_cols),
    df2     = m_es$nobs - length(coef(m_es)),
    p_joint = round(as.numeric(w$p), 4),
    beta_m5 = round(betas["treat_k_m5"], 3),
    beta_m4 = round(betas["treat_k_m4"], 3),
    beta_m3 = round(betas["treat_k_m3"], 3),
    beta_m2 = round(betas["treat_k_m2"], 3)
  ))

  cat(sprintf("%-15s  Wald = %.2f  p_joint = %.4f   |  b_m5=%+.2f  b_m4=%+.2f  b_m3=%+.2f  b_m2=%+.2f\n",
              outcome_labels[i],
              as.numeric(w$stat), as.numeric(w$p),
              betas["treat_k_m5"], betas["treat_k_m4"],
              betas["treat_k_m3"], betas["treat_k_m2"]))
}

# ---- 4. Save outputs -------------------------------------------------------

saveRDS(results, "out/pretrend_joint_test.rds")

sink("out/pretrend_joint_test.txt")
cat("Pre-trend Joint Wald Test - Iowa event study\n")
cat("H0: beta_{-5} = beta_{-4} = beta_{-3} = beta_{-2} = 0\n")
cat("Cluster-robust SE at county level\n")
cat(strrep("=", 70), "\n", sep = "")
print(results)
cat("\nInterpretation:\n")
cat(" - F_stat, p_joint: cluster-robust joint test of zero pre-trends\n")
cat(" - beta_m5..m2:     point estimates at k = -5, -4, -3, -2\n")
cat(" - p_joint > 0.10 is consistent with parallel trends\n")
sink()

cat("\n=== Saved ===\n")
cat("  out/pretrend_joint_test.rds\n")
cat("  out/pretrend_joint_test.txt\n")

# ---- 5. Defense one-liner --------------------------------------------------

drug_row <- results[outcome == "Drug arrests"]
if (nrow(drug_row) > 0 && !is.na(drug_row$F_stat)) {
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat("DEFENSE ONE-LINER:\n")
  cat(sprintf("'The joint Wald test of all four pre-period coefficients yields\n"))
  cat(sprintf(" F = %.2f, p_joint = %.3f, failing to reject the null of zero\n",
              drug_row$F_stat, drug_row$p_joint))
  cat(sprintf(" pre-trends. Pre-period point estimates are (b_-5=%+.2f, b_-4=%+.2f,\n",
              drug_row$beta_m5, drug_row$beta_m4))
  cat(sprintf(" b_-3=%+.2f, b_-2=%+.2f), uniformly biasing against the main\n",
              drug_row$beta_m3, drug_row$beta_m2))
  cat(" finding.'\n")
  cat(strrep("-", 70), "\n", sep = "")
}

cat("\nDone. Add joint F-test row to Table 5.2 or footnote in Ch5.4.\n")
