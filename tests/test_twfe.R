# =============================================================================
# tests/test_twfe.R   (testthat)
#
# Sanity checks for R/utils/twfe_no_fixest.R. The acceptance test is that
# the fallback estimator returns coefficients consistent with fixest when
# fixest is available.
# =============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(data.table)
})
source(here::here("R", "utils", "twfe_no_fixest.R"))

# -----------------------------------------------------------------------------
# Synthetic panel
# -----------------------------------------------------------------------------
make_panel <- function(n_cty = 50, n_yr = 8, n_border = 10,
                       true_effect = 0.1, seed = 123) {
  set.seed(seed)
  panel <- CJ(county = sprintf("c%03d", seq_len(n_cty)),
              year = seq(2015, 2014 + n_yr))
  panel[, border := as.integer(as.integer(sub("c", "", county)) <= n_border)]
  panel[, post := as.integer(year >= 2020)]
  panel[, treat := border * post]
  panel[, ln_pop := log(50000 + rnorm(.N, 0, 5000))]
  panel[, y := 5 + true_effect * treat + rnorm(.N, 0, 0.3)]
  panel
}

test_that("fit_twfe recovers a planted +0.10 effect", {
  panel <- make_panel(true_effect = 0.10)
  fit <- fit_twfe(
    y_col = "y", treat_col = "treat", controls = "ln_pop",
    fe = c("county", "year"), data = panel, cluster_col = "county"
  )
  expect_true(abs(fit$beta - 0.10) < 0.05)
  expect_equal(fit$n_obs, 400L)
  expect_equal(fit$n_cluster, 50L)
})

test_that("fit_twfe returns near-zero when planted effect is zero", {
  panel <- make_panel(true_effect = 0)
  fit <- fit_twfe(
    y_col = "y", treat_col = "treat", controls = "ln_pop",
    fe = c("county", "year"), data = panel, cluster_col = "county"
  )
  expect_true(abs(fit$beta) < 0.05)
  expect_true(fit$p > 0.1)
})

test_that("fit_event_study returns symmetric output around k = -1", {
  panel <- make_panel()
  panel[, year_rel := year - 2020]
  es <- fit_event_study(
    y_col = "y", treat_col = "border",
    time_var = "year_rel", ref_k = -1L,
    controls = "ln_pop", fe = c("county", "year"),
    data = panel, cluster_col = "county"
  )
  expect_equal(nrow(es), 8L)
  ref_row <- es[k == -1]
  expect_equal(ref_row$est, 0)
  expect_true(is.na(ref_row$pval))
})

test_that("fixest::feols agrees with the lm fallback (when fixest installed)", {
  skip_if_not_installed("fixest")
  panel <- make_panel(true_effect = 0.07)
  fb <- fit_twfe(
    y_col = "y", treat_col = "treat", controls = "ln_pop",
    fe = c("county", "year"), data = panel, cluster_col = "county"
  )
  m <- fixest::feols(y ~ treat + ln_pop | county + year, data = panel,
                     cluster = ~ county)
  fx_beta <- coef(m)["treat"]
  fx_se   <- sqrt(vcov(m)["treat", "treat"])
  expect_true(abs(fb$beta - fx_beta) < 1e-4)
  expect_true(abs(fb$se - fx_se) < 1e-3)
})
