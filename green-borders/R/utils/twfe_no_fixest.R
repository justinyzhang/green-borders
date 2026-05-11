# =============================================================================
# utils/twfe_no_fixest.R
#
# TWFE estimation without fixest, using base lm() + sandwich + lmtest.
#
# Two public functions:
#   fit_twfe()         — single TWFE regression with cluster-robust SE
#   fit_event_study()  — event study with k bin dummies × treatment
#
# Both return a tidy list mirroring the parts of the fixest API we use:
#   list(beta, se, p, ci_lo, ci_hi, n_obs, n_cluster, model, vcov)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
})

#' Fit a TWFE regression with cluster-robust SE
#'
#' @param y_col character. Outcome column name.
#' @param treat_col character. Treatment column name (the coefficient of interest).
#' @param controls character vector. Additional control columns.
#' @param fe character vector. Fixed-effect column names. Factorized inline.
#' @param data data.frame or data.table. Must contain all referenced columns.
#' @param cluster_col character. Column for cluster-robust SE.
#' @param weights numeric or NULL. Observation weights.
#'
#' @return list with elements beta, se, p, ci_lo, ci_hi, n_obs, n_cluster,
#'   model, vcov.
fit_twfe <- function(y_col, treat_col, controls = character(0),
                     fe = c("county_fips", "year"),
                     data, cluster_col = "county_fips",
                     weights = NULL) {
  stopifnot(y_col %in% names(data),
            treat_col %in% names(data),
            all(controls %in% names(data)),
            all(fe %in% names(data)),
            cluster_col %in% names(data))
  
  rhs_terms <- c(treat_col, controls,
                 sprintf("factor(%s)", fe))
  formula_str <- sprintf("%s ~ %s", y_col, paste(rhs_terms, collapse = " + "))
  f <- as.formula(formula_str)
  
  dt <- as.data.frame(data)
  if (is.null(weights)) {
    m <- lm(f, data = dt)
  } else {
    m <- lm(f, data = dt, weights = weights)
  }
  
  vc <- sandwich::vcovCL(m, cluster = dt[[cluster_col]], type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  
  list(
    beta      = unname(ct[treat_col, "Estimate"]),
    se        = unname(ct[treat_col, "Std. Error"]),
    p         = unname(ct[treat_col, "Pr(>|t|)"]),
    ci_lo     = unname(ci[treat_col, 1]),
    ci_hi     = unname(ci[treat_col, 2]),
    n_obs     = nobs(m),
    n_cluster = length(unique(dt[[cluster_col]])),
    model     = m,
    vcov      = vc
  )
}

#' Event study with k bin dummies x treatment
#'
#' @param y_col character. Outcome column name.
#' @param treat_col character. Treatment indicator (border / 0-1).
#' @param time_var character. Continuous event time, e.g. year_rel = year - 2020.
#' @param ref_k integer. Reference period (typically -1).
#' @param fe,controls,data,cluster_col,weights as in fit_twfe.
#'
#' @return data.table with columns k, est, se, ci_lo, ci_hi, pval. The
#'   reference row (k = ref_k) has all zeros and NA pval.
fit_event_study <- function(y_col, treat_col = "border",
                            time_var = "year_rel", ref_k = -1L,
                            controls = "ln_pop",
                            fe = c("county_fips", "year"),
                            data, cluster_col = "county_fips",
                            weights = NULL) {
  dt <- as.data.table(copy(data))
  ks <- sort(unique(dt[[time_var]]))
  for (k in ks) {
    if (k == ref_k) next
    cname <- sprintf("k_%s", ifelse(k < 0, paste0("m", abs(k)), as.character(k)))
    dt[, (cname) := as.integer(get(time_var) == k) * get(treat_col)]
  }
  k_cols <- grep("^k_", names(dt), value = TRUE)
  
  rhs_terms <- c(k_cols, controls, sprintf("factor(%s)", fe))
  formula_str <- sprintf("%s ~ %s", y_col, paste(rhs_terms, collapse = " + "))
  f <- as.formula(formula_str)
  
  df <- as.data.frame(dt)
  if (is.null(weights)) {
    m <- lm(f, data = df)
  } else {
    m <- lm(f, data = df, weights = weights)
  }
  
  vc <- sandwich::vcovCL(m, cluster = df[[cluster_col]], type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  
  rows <- data.table(
    k     = sapply(k_cols, function(s) {
      x <- sub("^k_", "", s)
      if (substr(x, 1, 1) == "m") -as.integer(sub("^m", "", x)) else as.integer(x)
    }),
    est   = ct[k_cols, "Estimate"],
    se    = ct[k_cols, "Std. Error"],
    pval  = ct[k_cols, "Pr(>|t|)"],
    ci_lo = ci[k_cols, 1],
    ci_hi = ci[k_cols, 2]
  )
  rows <- rbindlist(list(rows, data.table(k = ref_k, est = 0, se = 0,
                                          pval = NA_real_, ci_lo = 0, ci_hi = 0)))
  setorder(rows, k)
  rows[]
}
