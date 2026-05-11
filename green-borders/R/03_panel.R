# =============================================================================
# 03_panel.R   (4-state expansion: IA + WI + IN + MO)
#
# Build analytical panel:
#   - county_fips x state_name x year
#   - treatment status (IL_border, interior, MO_IL_border, MO_interior)
#   - NIBRS arrest counts (drug, owi, property, violent)
#   - ACS demographics
#   - DiD variables
#
# Outputs:
#   data/processed/panel_county_year.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

meta <- fread(here::here("data", "processed", "iowa_county_metadata.csv"))
# Ensure county_fips is 5-char zero-padded character
meta[, county_fips := sprintf("%05s", as.character(county_fips))]
cat("Metadata:", nrow(meta), "counties\n")

nibrs <- readRDS(here::here("data", "interim", "nibrs_county_year.rds"))
setDT(nibrs)
nibrs[, county_fips := as.character(county_fips)]
cat("NIBRS:   ", nrow(nibrs), "county-year cells\n")

acs <- readRDS(here::here("data", "interim", "acs_iowa.rds"))
setDT(acs)
acs[, county_fips := as.character(county_fips)]
cat("ACS:     ", nrow(acs), "county-year cells\n")

all_years <- 2015:2022
panel <- CJ(county_fips = meta$county_fips, year = all_years)
cat("Skeleton:", nrow(panel), "rows\n")

# Merge metadata
panel <- merge(panel,
               meta[, .(county_fips, name, state_name, state_fips,
                        border, mo_border, treatment, border_via,
                        dist_to_il_border_mi, land_area_sqmi)],
               by = "county_fips", all.x = TRUE)

# Merge NIBRS
panel <- merge(panel, nibrs[, .(county_fips, year, drug, owi, property, violent)],
               by = c("county_fips", "year"), all.x = TRUE)
for (col in c("owi", "drug", "property", "violent")) {
  if (!col %in% names(panel)) panel[, (col) := 0L]
  panel[is.na(get(col)), (col) := 0L]
}

# Merge ACS
panel <- merge(panel, acs, by = c("county_fips", "year"), all.x = TRUE)
n_missing_pop <- sum(is.na(panel$pop_total))
if (n_missing_pop > 0) cat("WARN:", n_missing_pop, "county-years missing pop\n")

# DiD variables: treat = IL-border (excluding Missouri) x post-2020
panel[, post := as.integer(year >= 2020)]
panel[, is_il_border := as.integer(border == 1 & state_name != "missouri")]
panel[, treat := is_il_border * post]
panel[, year_rel := year - 2020]

# Missouri-specific: 2023 placebo (NA in 2015-2022 window since panel ends 2022,
# but flag MO-IL border for separate analysis)
panel[, mo_post_2023 := as.integer(state_name == "missouri" & year >= 2023)]
panel[, mo_treat := as.integer(border == 1 & state_name == "missouri") * mo_post_2023]

# Rates and logs
for (col in c("owi", "drug", "property", "violent")) {
  rate_col <- paste0(col, "_rate")
  log_col  <- paste0("ln_", col)
  panel[, (rate_col) := get(col) / pop_total * 1e5]
  panel[, (log_col)  := log(get(col) + 1)]
}
panel[, ln_pop := log(pop_total)]
panel[, ln_inc := log(med_inc)]
panel[, pov_rate := pov_total / pop_total]

# Factor variables for FE
panel[, state_year := paste0(state_name, "_", year)]   # for state-by-year FE
panel[, state_name := factor(state_name)]

setorder(panel, county_fips, year)
dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(panel, here::here("data", "processed", "panel_county_year.rds"))

cat("\n=== Panel summary ===\n")
cat("Rows total:        ", nrow(panel), "\n")
cat("Unique counties:   ", length(unique(panel$county_fips)), "\n")
cat("Year range:        ", paste(range(panel$year), collapse = " - "), "\n")
cat("\n=== By treatment x state ===\n")
print(panel[, .N, by = .(state_name, treatment)])

cat("\n=== Treated obs (IL_border x post) ===\n")
cat("Total treat==1 obs:", sum(panel$treat), "\n")
cat("Distinct IL-border treated counties:",
    length(unique(panel[is_il_border == 1]$county_fips)), "\n")

cat("\n=== Mean count by treatment x post ===\n")
print(panel[, .(owi      = round(mean(owi, na.rm = TRUE), 1),
                drug     = round(mean(drug, na.rm = TRUE), 1),
                property = round(mean(property, na.rm = TRUE), 1),
                violent  = round(mean(violent, na.rm = TRUE), 1),
                pop      = round(mean(pop_total, na.rm = TRUE), 0),
                n        = .N),
            by = .(treatment, post)])

cat("\n=== Naive DiD (rate change, IL_border vs interior, excluding MO) ===\n")
diff_rates <- function(outcome_col) {
  rate_col <- paste0(outcome_col, "_rate")
  treated_set <- panel[treatment == "IL_border"]
  ctrl_set    <- panel[treatment == "interior"]
  tp <- mean(treated_set[post == 0][[rate_col]], na.rm = TRUE)
  tq <- mean(treated_set[post == 1][[rate_col]], na.rm = TRUE)
  cp <- mean(ctrl_set[post == 0][[rate_col]], na.rm = TRUE)
  cq <- mean(ctrl_set[post == 1][[rate_col]], na.rm = TRUE)
  did <- (tq - tp) - (cq - cp)
  data.table(outcome = outcome_col,
             treat_pre = round(tp, 1), treat_post = round(tq, 1),
             ctrl_pre = round(cp, 1), ctrl_post = round(cq, 1),
             naive_did = round(did, 1), pct_did = round(100 * did / tp, 1))
}
print(rbindlist(lapply(c("owi","drug","property","violent"), diff_rates)))

cat("\nWrote panel_county_year.rds\n")
