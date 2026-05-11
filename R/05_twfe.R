# =============================================================================
# 05_twfe.R   (4-state: IA + WI + IN; MO excluded from main DiD)
#
# Baseline TWFE with cluster-robust SE at county level (sandwich HC1).
#
# Specification (Sergey-recommended):
#   ln(Y_cst + 1) = alpha_c + lambda_st + tau * (Border_c x Post_t)
#                + gamma_1 * ln_pop_cst + gamma_2 * ln_inc_cst + u_cst
#
# where lambda_st is a STATE-BY-YEAR fixed effect (absorbs each state's
# specific time trajectory; allows IA / IN / WI to have their own panels).
#
# Outputs:
#   out/tab1_baseline_coefs.csv
#   out/tab1_baseline.tex
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(here)
})

panel <- readRDS(here::here("data", "processed", "panel_county_year.rds"))
setDT(panel)

# Main DiD sample: 3 prohibition states (excluding Missouri)
panel_main <- panel[treatment %in% c("IL_border", "interior")]
cat("Main DiD sample:", nrow(panel_main), "rows,",
    length(unique(panel_main$county_fips)), "counties\n")
cat("Treated counties:", length(unique(panel_main[is_il_border == 1]$county_fips)), "\n")
cat("State distribution of treated:\n")
print(panel_main[is_il_border == 1 & year == 2020, .N, by = state_name])

fit_twfe <- function(y_col, dt = panel_main, with_state_year_fe = TRUE) {
  if (with_state_year_fe) {
    f <- as.formula(sprintf(
      "%s ~ treat + ln_pop + ln_inc + factor(county_fips) + factor(state_year)",
      y_col
    ))
  } else {
    f <- as.formula(sprintf(
      "%s ~ treat + ln_pop + ln_inc + factor(county_fips) + factor(year)",
      y_col
    ))
  }
  m <- lm(f, data = dt)
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  list(model = m, ct = ct, ci = ci)
}

cat("\n=== Spec A: county FE + state-by-year FE (preferred) ===\n")
fit_drug_A     <- fit_twfe("ln_drug",     with_state_year_fe = TRUE)
fit_owi_A      <- fit_twfe("ln_owi",      with_state_year_fe = TRUE)
fit_property_A <- fit_twfe("ln_property", with_state_year_fe = TRUE)
fit_violent_A  <- fit_twfe("ln_violent",  with_state_year_fe = TRUE)

cat("\n=== Spec B: county FE + year FE (no state interaction) ===\n")
fit_drug_B     <- fit_twfe("ln_drug",     with_state_year_fe = FALSE)
fit_owi_B      <- fit_twfe("ln_owi",      with_state_year_fe = FALSE)
fit_property_B <- fit_twfe("ln_property", with_state_year_fe = FALSE)
fit_violent_B  <- fit_twfe("ln_violent",  with_state_year_fe = FALSE)

make_row <- function(fit, label, spec) {
  b   <- fit$ct["treat", "Estimate"]
  se  <- fit$ct["treat", "Std. Error"]
  pv  <- fit$ct["treat", "Pr(>|t|)"]
  ci  <- fit$ci["treat", ]
  data.table(
    outcome = label, spec = spec,
    beta    = round(b, 4),
    se      = round(se, 4),
    p_value = round(pv, 4),
    ci_lo   = round(ci[1], 4),
    ci_hi   = round(ci[2], 4),
    n_obs   = nobs(fit$model)
  )
}

coef_summary <- rbindlist(list(
  make_row(fit_drug_A,     "Drug",     "A_state_year_FE"),
  make_row(fit_owi_A,      "OWI",      "A_state_year_FE"),
  make_row(fit_property_A, "Property", "A_state_year_FE"),
  make_row(fit_violent_A,  "Violent",  "A_state_year_FE"),
  make_row(fit_drug_B,     "Drug",     "B_year_FE"),
  make_row(fit_owi_B,      "OWI",      "B_year_FE"),
  make_row(fit_property_B, "Property", "B_year_FE"),
  make_row(fit_violent_B,  "Violent",  "B_year_FE")
))

cat("\n=== Baseline TWFE summary ===\n")
print(coef_summary)

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)
fwrite(coef_summary, here::here("out", "tab1_baseline_coefs.csv"))
cat("\nWrote out/tab1_baseline_coefs.csv\n")

# -----------------------------------------------------------------------------
# LaTeX (4-column: 4 outcomes, spec A only)
# -----------------------------------------------------------------------------
stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("$^{***}$")
  if (p < 0.05) return("$^{**}$")
  if (p < 0.10) return("$^{*}$")
  ""
}
fmt_b  <- function(b, p) sprintf("%.3f%s", b, stars(p))
fmt_se <- function(se) sprintf("(%.3f)", se)

A <- coef_summary[spec == "A_state_year_FE"]

tex_lines <- c(
  "\\begin{table}[!ht]\\centering",
  "\\caption{Cross-border spillovers of Illinois cannabis legalization on IL-bordering counties in three prohibition states (IA, IN, WI)}",
  "\\label{tab:baseline}",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "& Drug & OWI & Property & Violent \\\\",
  "& (1) & (2) & (3) & (4) \\\\",
  "\\midrule",
  paste0("Border $\\times$ Post & ",
         fmt_b(A[outcome == "Drug"]$beta,     A[outcome == "Drug"]$p_value),     " & ",
         fmt_b(A[outcome == "OWI"]$beta,      A[outcome == "OWI"]$p_value),      " & ",
         fmt_b(A[outcome == "Property"]$beta, A[outcome == "Property"]$p_value), " & ",
         fmt_b(A[outcome == "Violent"]$beta,  A[outcome == "Violent"]$p_value),  " \\\\"),
  paste0(" & ",
         fmt_se(A[outcome == "Drug"]$se),     " & ",
         fmt_se(A[outcome == "OWI"]$se),      " & ",
         fmt_se(A[outcome == "Property"]$se), " & ",
         fmt_se(A[outcome == "Violent"]$se),  " \\\\"),
  "\\midrule",
  "County FE      & Yes & Yes & Yes & Yes \\\\",
  "State-by-Year FE & Yes & Yes & Yes & Yes \\\\",
  "Controls       & Yes & Yes & Yes & Yes \\\\",
  paste0("Observations & ",
         A[outcome == "Drug"]$n_obs,     " & ",
         A[outcome == "OWI"]$n_obs,      " & ",
         A[outcome == "Property"]$n_obs, " & ",
         A[outcome == "Violent"]$n_obs,  " \\\\"),
  paste0("Treated counties & 27 & 27 & 27 & 27 \\\\"),
  paste0("Control counties & 236 & 236 & 236 & 236 \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{0.95\\textwidth}\\footnotesize",
  "\\textit{Notes}: Outcome variables are $\\log(\\mathrm{arrests} + 1)$. The sample comprises 263 counties in Iowa, Indiana, and Wisconsin observed from 2015 to 2022. Treated counties (27) are those directly bordering Illinois; control counties (236) are interior counties in the same three states. All specifications include county fixed effects, state-by-year fixed effects, log population, and log median household income. Missouri is excluded from the main panel and analyzed separately as a placebo (Table B). Cluster-robust standard errors at the county level (HC1) in parentheses. $^{*}$ $p < 0.10$, $^{**}$ $p < 0.05$, $^{***}$ $p < 0.01$.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tex_lines, here::here("out", "tab1_baseline.tex"))
cat("Wrote out/tab1_baseline.tex\n")
