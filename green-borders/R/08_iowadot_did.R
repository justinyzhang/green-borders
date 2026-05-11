# =============================================================================
# 08_iowadot_did.R   (DiD on Iowa DOT outcomes, Iowa-only, 2015-2022)
#
# Re-runs the same Iowa-only DiD setup using THREE cleaner Iowa DOT outcomes
# instead of NIBRS arrests:
#   1. OWI revocations (direct legal outcome)
#   2. Speeding convictions (enforcement intensity proxy)
#   3. Yearly fatalities (downstream harm)
#
# These are state government registry data with NO NIBRS reporting volatility.
#
# Same DiD specification:
#   ln(Y_ct + 1) = alpha_c + lambda_t + tau * (Border_c x Post_t)
#                + ln_pop_ct + ln_inc_ct + u_ct
# Cluster-robust SE at county level.
#
# Output:
#   out/tab3_iowadot_baseline.csv
#   out/fig4_iowadot_eventstudy.pdf
#   out/tab3_iowadot_robust.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(here)
})

# -----------------------------------------------------------------------------
# 1. Build Iowa DOT panel using existing metadata + ACS
# -----------------------------------------------------------------------------
meta <- fread(here::here("data", "processed", "iowa_county_metadata.csv"))
meta[, county_fips := sprintf("%05s", as.character(county_fips))]
# If meta is from 4-state, restrict to Iowa
meta_ia <- meta[state_name == "iowa" | substr(county_fips, 1, 2) == "19"]

iadot <- readRDS(here::here("data", "interim", "iowadot_county_year.rds"))
setDT(iadot)
iadot[, county_fips := as.character(county_fips)]

acs <- readRDS(here::here("data", "interim", "acs_iowa.rds"))
setDT(acs)
acs[, county_fips := as.character(county_fips)]
acs_ia <- acs[substr(county_fips, 1, 2) == "19"]

all_counties <- meta_ia$county_fips
panel <- CJ(county_fips = all_counties, year = 2015:2022)

panel <- merge(panel,
               meta_ia[, .(county_fips, name, border, mo_border,
                            dist_to_il_border_mi, land_area_sqmi)],
               by = "county_fips", all.x = TRUE)
panel <- merge(panel, iadot, by = c("county_fips", "year"), all.x = TRUE)
panel <- merge(panel, acs_ia, by = c("county_fips", "year"), all.x = TRUE)

# DiD vars
panel[, post := as.integer(year >= 2020)]
panel[, treat := border * post]
panel[, year_rel := year - 2020]
panel[, ln_pop := log(pop_total)]
panel[, ln_inc := log(med_inc)]

# Log transforms for each available outcome
outcome_cols <- intersect(c("owi_revocation", "speeding_conviction", "fatality"),
                           names(panel))
for (col in outcome_cols) {
  rate_col <- paste0(col, "_rate")
  log_col <- paste0("ln_", col)
  panel[, (rate_col) := get(col) / pop_total * 1e5]
  panel[, (log_col) := log(pmax(get(col), 0) + 1)]
}

cat("=== Iowa DOT panel ===\n")
cat("Rows:", nrow(panel), "\n")
cat("Counties:", uniqueN(panel$county_fips), "\n")
cat("Outcomes available:", paste(outcome_cols, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# 2. Naive DiD (no FE, raw rates) - direction check
# -----------------------------------------------------------------------------
cat("=== Naive DiD (per-100k rates) ===\n")
diff_rates <- function(outcome_col) {
  rate_col <- paste0(outcome_col, "_rate")
  if (!rate_col %in% names(panel)) return(NULL)
  tp <- mean(panel[border == 1 & post == 0][[rate_col]], na.rm = TRUE)
  tq <- mean(panel[border == 1 & post == 1][[rate_col]], na.rm = TRUE)
  cp <- mean(panel[border == 0 & post == 0][[rate_col]], na.rm = TRUE)
  cq <- mean(panel[border == 0 & post == 1][[rate_col]], na.rm = TRUE)
  did <- (tq - tp) - (cq - cp)
  data.table(outcome = outcome_col,
             treat_pre = round(tp, 1), treat_post = round(tq, 1),
             ctrl_pre = round(cp, 1), ctrl_post = round(cq, 1),
             naive_did = round(did, 1),
             pct_did = round(100 * did / tp, 1))
}
naive <- rbindlist(lapply(outcome_cols, diff_rates))
print(naive)

# -----------------------------------------------------------------------------
# 3. TWFE with cluster SE for each outcome
# -----------------------------------------------------------------------------
fit_twfe <- function(y_col, dt = panel) {
  if (!y_col %in% names(dt)) return(NULL)
  f <- as.formula(sprintf(
    "%s ~ treat + ln_pop + ln_inc + factor(county_fips) + factor(year)",
    y_col
  ))
  m <- lm(f, data = dt[!is.na(get(y_col))])
  vc <- sandwich::vcovCL(m,
                          cluster = dt[!is.na(get(y_col))]$county_fips,
                          type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  list(beta = ct["treat","Estimate"],
       se = ct["treat","Std. Error"],
       p = ct["treat","Pr(>|t|)"],
       ci_lo = ci["treat",1], ci_hi = ci["treat",2],
       n = nobs(m))
}

cat("\n=== TWFE baseline ===\n")
twfe_rows <- list()
for (col in outcome_cols) {
  ln_col <- paste0("ln_", col)
  fit <- fit_twfe(ln_col)
  if (!is.null(fit)) {
    twfe_rows[[length(twfe_rows) + 1]] <- data.table(
      outcome = col, beta = round(fit$beta, 4), se = round(fit$se, 4),
      p_value = round(fit$p, 4),
      ci_lo = round(fit$ci_lo, 4), ci_hi = round(fit$ci_hi, 4),
      n_obs = fit$n
    )
  }
}
twfe_summary <- rbindlist(twfe_rows)
print(twfe_summary)

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)
fwrite(twfe_summary, here::here("out", "tab3_iowadot_baseline.csv"))
cat("Wrote out/tab3_iowadot_baseline.csv\n")

# -----------------------------------------------------------------------------
# 4. Event study on the strongest outcome (loop over all 3)
# -----------------------------------------------------------------------------
event_years <- sort(unique(panel$year_rel))
ref_k <- -1

fit_event <- function(y_col) {
  if (!y_col %in% names(panel)) return(NULL)
  dt <- copy(panel[!is.na(get(y_col))])
  for (k in event_years) {
    if (k == ref_k) next
    cname <- sprintf("k_%s", ifelse(k < 0, paste0("m", abs(k)), as.character(k)))
    dt[, (cname) := as.integer(year_rel == k) * border]
  }
  k_cols <- grep("^k_", names(dt), value = TRUE)
  f <- as.formula(sprintf(
    "%s ~ %s + ln_pop + factor(county_fips) + factor(year)",
    y_col, paste(k_cols, collapse = " + ")
  ))
  m <- lm(f, data = dt)
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  rows <- data.table(
    k = sapply(k_cols, function(s) {
      x <- sub("^k_", "", s)
      if (substr(x, 1, 1) == "m") -as.integer(sub("^m", "", x)) else as.integer(x)
    }),
    est = ct[k_cols, "Estimate"],
    se = ct[k_cols, "Std. Error"],
    ci_lo = ci[k_cols, 1], ci_hi = ci[k_cols, 2]
  )
  rbindlist(list(rows, data.table(k = ref_k, est = 0, se = 0,
                                   ci_lo = 0, ci_hi = 0)))[order(k)]
}

cat("\n=== Event studies ===\n")
all_es <- list()
for (col in outcome_cols) {
  cat("\nOutcome:", col, "\n")
  es <- fit_event(paste0("ln_", col))
  if (!is.null(es)) {
    print(es)
    all_es[[col]] <- data.table(outcome = col, es)
  }
}
if (length(all_es) > 0) {
  fwrite(rbindlist(all_es), here::here("out", "fig4_iowadot_eventstudy_coefs.csv"))
  
  combined_es <- rbindlist(all_es)
  combined_es[, outcome_label := fcase(
    outcome == "owi_revocation",       "OWI revocations",
    outcome == "speeding_conviction",  "Speeding convictions (enforcement intensity)",
    outcome == "fatality",             "Traffic fatalities"
  )]
  
  g <- ggplot(combined_es, aes(k, est)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15,
                  color = "#1f78b4", linewidth = 0.6) +
    geom_point(size = 2, color = "#1f78b4") +
    facet_wrap(~ outcome_label, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = event_years) +
    labs(x = "Years relative to January 2020 IL legalization",
         y = "Coefficient (95% CI)",
         title = "Iowa DOT outcomes: event study (border vs interior counties)",
         subtitle = "Iowa DOT registry data; reference year 2019") +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
  
  ggsave(here::here("out", "fig4_iowadot_eventstudy.pdf"),
         g, width = 7, height = 8)
  cat("\nWrote out/fig4_iowadot_eventstudy.pdf\n")
}

# -----------------------------------------------------------------------------
# 5. Robustness: drop Scott, drop 2020 COVID
# -----------------------------------------------------------------------------
cat("\n=== Robustness (5 specs each outcome) ===\n")
specs <- list(
  baseline    = panel,
  drop_2020   = panel[year != 2020],
  drop_scott  = panel[county_fips != "19163"],
  drop_mo     = panel[mo_border == 0]
)
robust_rows <- list()
for (col in outcome_cols) {
  ln_col <- paste0("ln_", col)
  for (sp_nm in names(specs)) {
    fit <- fit_twfe(ln_col, specs[[sp_nm]])
    if (!is.null(fit)) {
      robust_rows[[length(robust_rows) + 1]] <- data.table(
        outcome = col, spec = sp_nm,
        beta = round(fit$beta, 4), se = round(fit$se, 4),
        p_value = round(fit$p, 4), n = fit$n
      )
    }
  }
}
robust <- rbindlist(robust_rows)
print(robust)
fwrite(robust, here::here("out", "tab3_iowadot_robust.csv"))
cat("Wrote out/tab3_iowadot_robust.csv\n")
