# =============================================================================
# 07_robust.R   (4-state robustness on DRUG outcome)
#
# Six robustness specifications:
#   (1) Baseline (county FE + state-by-year FE)
#   (2) Drop 2020 (COVID donut)
#   (3) Drop Lake County, IN (Chicago suburb spillover from IL)
#   (4) Drop Scott County, IA (Quad Cities outlier)
#   (5) Drop all border-via-Mississippi-River (keep WI + IN land-border only)
#   (6) Population-weighted
#   (7) Within-state DiD: each state run separately (IA-only, IN-only, WI-only)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(here)
})

panel <- readRDS(here::here("data", "processed", "panel_county_year.rds"))
setDT(panel)
panel_main <- panel[treatment %in% c("IL_border", "interior")]

fit_twfe <- function(dt, y_col = "ln_drug", weighted = FALSE) {
  formula_text <- sprintf(
    "%s ~ treat + ln_pop + ln_inc + factor(county_fips) + factor(state_year)",
    y_col
  )
  if (weighted) {
    m <- lm(as.formula(formula_text), data = dt, weights = dt$pop_total)
  } else {
    m <- lm(as.formula(formula_text), data = dt)
  }
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  list(beta = ct["treat", "Estimate"], se = ct["treat", "Std. Error"],
       p = ct["treat", "Pr(>|t|)"],
       ci_lo = ci["treat", 1], ci_hi = ci["treat", 2],
       n = nobs(m), n_cluster = length(unique(dt$county_fips)))
}

# Within-state DiD (no state-by-year FE since single state)
fit_within_state <- function(dt, y_col = "ln_drug") {
  m <- lm(as.formula(sprintf(
    "%s ~ treat + ln_pop + ln_inc + factor(county_fips) + factor(year)",
    y_col)), data = dt)
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  list(beta = ct["treat", "Estimate"], se = ct["treat", "Std. Error"],
       p = ct["treat", "Pr(>|t|)"],
       ci_lo = ci["treat", 1], ci_hi = ci["treat", 2],
       n = nobs(m), n_cluster = length(unique(dt$county_fips)))
}

# Lake County IN FIPS = 18089; Scott County IA = 19163
specs <- list(
  baseline    = list(dt = panel_main,                              fitter = fit_twfe),
  drop_2020   = list(dt = panel_main[year != 2020],                fitter = fit_twfe),
  drop_lake   = list(dt = panel_main[county_fips != "18089"],      fitter = fit_twfe),
  drop_scott  = list(dt = panel_main[county_fips != "19163"],      fitter = fit_twfe),
  pop_weight  = list(dt = panel_main,                              fitter = function(d) fit_twfe(d, weighted = TRUE)),
  IA_only     = list(dt = panel_main[state_name == "iowa"],        fitter = fit_within_state),
  IN_only     = list(dt = panel_main[state_name == "indiana"],     fitter = fit_within_state),
  WI_only     = list(dt = panel_main[state_name == "wisconsin"],   fitter = fit_within_state)
)

results <- rbindlist(lapply(names(specs), function(nm) {
  out <- specs[[nm]]$fitter(specs[[nm]]$dt)
  data.table(
    spec      = nm,
    beta      = round(out$beta, 4),
    se        = round(out$se, 4),
    p_value   = round(out$p, 4),
    ci_lo     = round(out$ci_lo, 4),
    ci_hi     = round(out$ci_hi, 4),
    n_obs     = out$n,
    n_cluster = out$n_cluster
  )
}))

cat("\n=== Robustness: drug arrests (4-state, IL-border DiD) ===\n")
print(results)

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)
fwrite(results, here::here("out", "tab2_robust_coefs.csv"))

# LaTeX
stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("$^{***}$")
  if (p < 0.05) return("$^{**}$")
  if (p < 0.10) return("$^{*}$")
  ""
}
fmt_b  <- function(b, p) sprintf("%.3f%s", b, stars(p))
fmt_se <- function(se) sprintf("(%.3f)", se)

main_specs <- results[spec %in% c("baseline", "drop_2020", "drop_lake",
                                   "drop_scott", "pop_weight")]
tex_lines <- c(
  "\\begin{table}[!ht]\\centering",
  "\\caption{Robustness: drug arrests, alternative samples and weighting}",
  "\\label{tab:robust}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "& Baseline & Drop 2020 & Drop Lake (IN) & Drop Scott (IA) & Pop.-wt. \\\\",
  "& (1) & (2) & (3) & (4) & (5) \\\\",
  "\\midrule",
  paste0("Border $\\times$ Post & ",
         paste(sapply(seq_len(nrow(main_specs)), function(i)
                fmt_b(main_specs$beta[i], main_specs$p_value[i])), collapse = " & "),
         " \\\\"),
  paste0(" & ",
         paste(sapply(seq_len(nrow(main_specs)), function(i)
                fmt_se(main_specs$se[i])), collapse = " & "),
         " \\\\"),
  "\\midrule",
  paste0("Observations & ",
         paste(main_specs$n_obs, collapse = " & "),
         " \\\\"),
  paste0("Counties & ",
         paste(main_specs$n_cluster, collapse = " & "),
         " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.95\\textwidth}\\footnotesize",
  "\\textit{Notes}: Outcome is $\\log(\\mathrm{drug\\ arrests} + 1)$. Sample of IA, IN, and WI; Missouri excluded (analyzed separately). Column 1 is the baseline with county and state-by-year fixed effects. Column 2 drops 2020 (COVID donut). Column 3 drops Lake County, IN (Chicago suburb). Column 4 drops Scott County, IA (Quad Cities). Column 5 weights by county population. Cluster-robust SE at the county level. $^{*}$ $p<0.10$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tex_lines, here::here("out", "tab2_robust.tex"))
cat("Wrote out/tab2_robust.tex\n")
