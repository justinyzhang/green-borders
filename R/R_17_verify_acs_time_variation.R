# ============================================================================
# R/17_verify_acs_time_variation.R
# ----------------------------------------------------------------------------
# Verify whether ACS population and income controls in the Iowa panel
# are county-year time-varying (good) or county-only time-invariant (bad).
#
# Reviewer v9 concern: if ln_pop and ln_inc are only county-level (single
# ACS 2016-2020 cross-section), then county fixed effects ABSORB them
# completely, and they cannot serve as controls — meaning the spec is wrong.
#
# What we need: ln_pop_ct and ln_inc_ct varying within county across years.
#
# Inputs:  data/processed/panel_county_year.rds
# Outputs: console only (no need to save)
# ============================================================================

library(data.table)

panel <- readRDS("data/processed/panel_county_year.rds")
panel <- as.data.table(panel)

state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
ia_panel <- panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia")]

cat(sprintf("Iowa panel: %d county-year obs\n", nrow(ia_panel)))
cat(sprintf("Counties: %d, Years: %s\n",
            length(unique(ia_panel$county_fips)),
            paste(range(ia_panel$year), collapse = " to ")))

# Check column names to see what variables exist
cat("\n=== Column names in Iowa panel ===\n")
print(names(ia_panel))

# Test 1: Is ln_pop varying within county?
cat("\n=== TEST 1: Within-county variation in ln_pop ===\n")
var_test_pop <- ia_panel[, .(
  n_unique_pop = uniqueN(ln_pop),
  range_pop = max(ln_pop) - min(ln_pop)
), by = county_fips]

cat(sprintf("Counties with n_unique_pop = 1 (TIME-INVARIANT, BAD): %d\n",
            sum(var_test_pop$n_unique_pop == 1)))
cat(sprintf("Counties with n_unique_pop > 1 (varying, OK):         %d\n",
            sum(var_test_pop$n_unique_pop > 1)))
cat(sprintf("Mean log(pop) range within county:                    %.4f\n",
            mean(var_test_pop$range_pop, na.rm = TRUE)))

# Test 2: Is ln_inc varying within county?
cat("\n=== TEST 2: Within-county variation in ln_inc ===\n")
if ("ln_inc" %in% names(ia_panel)) {
  var_test_inc <- ia_panel[, .(
    n_unique_inc = uniqueN(ln_inc),
    range_inc = max(ln_inc) - min(ln_inc)
  ), by = county_fips]
  
  cat(sprintf("Counties with n_unique_inc = 1 (TIME-INVARIANT, BAD): %d\n",
              sum(var_test_inc$n_unique_inc == 1)))
  cat(sprintf("Counties with n_unique_inc > 1 (varying, OK):         %d\n",
              sum(var_test_inc$n_unique_inc > 1)))
  cat(sprintf("Mean log(inc) range within county:                    %.4f\n",
              mean(var_test_inc$range_inc, na.rm = TRUE)))
} else {
  cat("ln_inc column not found — check actual column name\n")
}

# Test 3: Show a sample county's ln_pop trajectory
cat("\n=== TEST 3: Sample county (Scott County, FIPS 19163) trajectory ===\n")
sample_county <- ia_panel[county_fips == "19163"]
if (nrow(sample_county) > 0) {
  print(sample_county[, .(year, ln_pop, ln_inc)])
} else {
  cat("Scott County not found — try another FIPS\n")
  print(ia_panel[1:8, .(county_fips, year, ln_pop)])
}

# ============================================================================
# VERDICT
# ============================================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("VERDICT\n")
cat(strrep("=", 70), "\n", sep = "")

pop_invariant_share <- sum(var_test_pop$n_unique_pop == 1) / nrow(var_test_pop)

if (pop_invariant_share > 0.5) {
  cat("⚠️  ln_pop appears TIME-INVARIANT (single ACS estimate per county).\n")
  cat("    In county FE specifications, ln_pop is FULLY ABSORBED.\n")
  cat("    The current thesis equation log(Pop_ct) is INCORRECT.\n\n")
  cat("    REQUIRED FIX in v13 thesis:\n")
  cat("    Option A: Re-extract data using annual ACS 5-year rolling estimates\n")
  cat("              (2015 panel = 2011-2015 ACS, 2016 = 2012-2016, etc.)\n")
  cat("    Option B: Remove ln_pop and ln_inc from controls, since county FE\n")
  cat("              absorbs them. Re-estimate baseline DiD without controls\n")
  cat("              and verify the headline β = 0.422 still holds (it should,\n")
  cat("              since FE absorbs time-invariant heterogeneity already).\n")
  cat("    Option C: Update Data 3.2 to clarify controls are time-invariant\n")
  cat("              county-level summaries that DROP OUT of county-FE specs.\n")
} else {
  cat("✓  ln_pop is COUNTY-YEAR TIME-VARYING.\n")
  cat("   Current thesis equation log(Pop_ct) is CORRECT.\n\n")
  cat("   Optional v13 addition: clarify in Data 3.2 that controls come from\n")
  cat("   ANNUAL ACS five-year rolling estimates (each year's 5-year window).\n")
}

cat("\nPaste this output to your Claude session.\n")
