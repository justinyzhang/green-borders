# =============================================================================
# 01_nibrs_load.R   (4-state expansion: IA + WI + IN + MO, no arrow, no fixest)
#
# Chapter 5 design (Sergey-recommended):
#   - Treated: IL-border counties in IA + WI + IN (and MO pre-2023)
#   - Control: same-state non-IL-border counties
#   - Shock:   Illinois January 2020 recreational legalization
#   - Outcome: NIBRS arrests (drug, OWI, property, violent)
#
# This script reads Kaplan V9 per-year .rds files, filters to 4 prohibition
# states, joins with LEAIC for county FIPS, and aggregates to county-year-
# offense_group counts. Missouri is included throughout; MO 2023 placebo
# logic is applied later in 03_panel.R / 07_robust.R.
#
# Inputs:
#   data/raw/nibrs/nibrs_{offense,administrative,group_b_arrest_report}_segment_{2015..2022}.rds
#   data/raw/leaic/35158-0001-Data.rda
#
# Output:
#   data/interim/nibrs_county_year.rds   (n_counties x 8 years long-format)
# =============================================================================

library(data.table)
library(here)

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Configuration: 4 prohibition states + IL for placebo / context
# -----------------------------------------------------------------------------
target_state_names <- c("iowa", "indiana", "wisconsin", "missouri")
target_state_abbs  <- c("ia", "in", "wi", "mo")
# Note: IL itself is excluded from the panel (it's the legalizing state)

DRUG_CODES <- c(
  "drug/narcotic offenses - drug/narcotic violations",
  "drug/narcotic offenses - drug equipment violations"
)
PROPERTY_CODES <- c(
  "burglary/breaking and entering",
  "larceny/theft offenses - all other larceny",
  "larceny/theft offenses - pocket-picking",
  "larceny/theft offenses - purse-snatching",
  "larceny/theft offenses - shoplifting",
  "larceny/theft offenses - theft from building",
  "larceny/theft offenses - theft from coin-operated machine or device",
  "larceny/theft offenses - theft from motor vehicle",
  "larceny/theft offenses - theft of motor vehicle parts/accessories",
  "motor vehicle theft",
  "destruction/damage/vandalism of property"
)
VIOLENT_CODES <- c(
  "assault offenses - aggravated assault",
  "assault offenses - simple assault",
  "assault offenses - intimidation",
  "robbery",
  "murder/nonnegligent manslaughter",
  "negligent manslaughter",
  "sex offenses - rape",
  "sex offenses - sodomy",
  "sex offenses - sexual assault with an object",
  "sex offenses - fondling (incident liberties/child molest)",
  "kidnapping/abduction"
)
OWI_CODES <- c(
  "driving under the influence",
  "driving under the influence of alcohol/drugs",
  "dui",
  "owi",
  "operating while intoxicated"
)
ALL_NONOWI <- c(DRUG_CODES, PROPERTY_CODES, VIOLENT_CODES)

target_years <- 2015:2022

# -----------------------------------------------------------------------------
# 2. LEAIC crosswalk: ORI -> 5-digit county FIPS for 4 target states
# -----------------------------------------------------------------------------
cat("[Setup] Loading LEAIC crosswalk ...\n")
leaic_path <- here::here("data", "raw", "leaic", "35158-0001-Data.rda")
if (!file.exists(leaic_path)) {
  cand <- list.files(here::here("data", "raw", "leaic"),
                     pattern = "\\.rda$", full.names = TRUE)
  if (length(cand) == 0) stop("Missing LEAIC at ", leaic_path)
  leaic_path <- cand[1]
}
leaic_env <- new.env()
load(leaic_path, envir = leaic_env)
leaic_obj <- ls(leaic_env)[1]
cat("  LEAIC object name:", leaic_obj, "\n")
leaic <- as.data.table(get(leaic_obj, envir = leaic_env))
rm(leaic_env)

for (col in c("STATENAME", "ORI9", "FIPS_ST", "FIPS_COUNTY")) {
  if (col %in% names(leaic)) leaic[[col]] <- as.character(leaic[[col]])
}

# 4-state subset
target_fips_codes <- c("18", "19", "29", "55")  # IN, IA, MO, WI
state_mask <- tolower(leaic$STATENAME) %in% target_state_names |
              leaic$FIPS_ST %in% target_fips_codes

leaic_4 <- leaic[state_mask & leaic$ORI9 != "-1" &
                 !is.na(leaic$ORI9) & nchar(leaic$ORI9) >= 7,
                 .(ori = ORI9,
                   state_name = tolower(STATENAME),
                   FIPS_ST, FIPS_COUNTY)]
leaic_4[, st_clean := gsub("\\D", "", FIPS_ST)]
leaic_4[, ct_clean := gsub("\\D", "", FIPS_COUNTY)]
leaic_4 <- leaic_4[nchar(st_clean) > 0 & nchar(ct_clean) > 0]
leaic_4[, county_fips := sprintf("%02d%03d",
                                  as.integer(st_clean),
                                  as.integer(ct_clean))]
leaic_4 <- unique(leaic_4[, .(ori, state_name, county_fips)], by = "ori")
cat("  ORIs across 4 states:", nrow(leaic_4), "\n")
print(leaic_4[, .N, by = state_name])

# -----------------------------------------------------------------------------
# 3. File discovery
# -----------------------------------------------------------------------------
nibrs_dir <- here::here("data", "raw", "nibrs")
all_files <- list.files(nibrs_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(all_files) == 0) stop("No .rds files in ", nibrs_dir)
cat("\nFound", length(all_files), ".rds files in nibrs/\n")

# -----------------------------------------------------------------------------
# 4. Per-year loop
# -----------------------------------------------------------------------------
collected <- list()

for (yr in target_years) {
  cat(sprintf("\n[%d] -----------------------\n", yr))
  
  # ---- offense segment ----
  off_path <- all_files[grepl(
    sprintf("offense_segment_%d", yr), basename(all_files)
  ) & !grepl("group_b", basename(all_files))]
  
  if (length(off_path) >= 1) {
    off <- readRDS(off_path[1])
    cat(sprintf("  offense file: %s (%s rows)\n",
                basename(off_path[1]),
                format(nrow(off), big.mark = ",")))
    
    if (is.factor(off$state))            off$state            <- as.character(off$state)
    if (is.factor(off$state_abb))        off$state_abb        <- as.character(off$state_abb)
    if (is.factor(off$ucr_offense_code)) off$ucr_offense_code <- as.character(off$ucr_offense_code)
    if (is.factor(off$ori))              off$ori              <- as.character(off$ori)
    
    target_state_mask <- tolower(off$state) %in% target_state_names |
                         tolower(off$state_abb) %in% target_state_abbs
    code_mask <- off$ucr_offense_code %in% ALL_NONOWI
    keep_mask <- target_state_mask & code_mask
    
    cat(sprintf("  offense: %s in 4 states, %s target code, %s both\n",
                format(sum(target_state_mask), big.mark = ","),
                format(sum(code_mask), big.mark = ","),
                format(sum(keep_mask), big.mark = ",")))
    
    if (sum(keep_mask) > 0) {
      collected[[length(collected) + 1]] <- data.table(
        ori = off$ori[keep_mask],
        offense_code = off$ucr_offense_code[keep_mask],
        year_inc = yr,
        source = "offense"
      )
    }
    rm(off); gc(verbose = FALSE)
  } else {
    cat("  offense: no file for", yr, "\n")
  }
  
  # ---- group_b arrest report (OWI) ----
  gb_path <- all_files[grepl(
    sprintf("group_b.*%d", yr), basename(all_files)
  )]
  
  if (length(gb_path) >= 1) {
    gb <- readRDS(gb_path[1])
    cat(sprintf("  group_b file: %s (%s rows)\n",
                basename(gb_path[1]),
                format(nrow(gb), big.mark = ",")))
    
    code_col_gb <- intersect(c("ucr_arrest_offense_code", "arrest_offense_code",
                                "ucr_offense_code", "offense_code"),
                              names(gb))[1]
    if (is.na(code_col_gb)) {
      cat("  WARN: no offense code column in group_b\n")
    } else {
      cat("  group_b code col:", code_col_gb, "\n")
      if (is.factor(gb$state))           gb$state           <- as.character(gb$state)
      if (is.factor(gb$state_abb))       gb$state_abb       <- as.character(gb$state_abb)
      if (is.factor(gb$ori))             gb$ori             <- as.character(gb$ori)
      if (is.factor(gb[[code_col_gb]]))  gb[[code_col_gb]]  <- as.character(gb[[code_col_gb]])
      
      target_state_mask_gb <- tolower(gb$state) %in% target_state_names |
                              tolower(gb$state_abb) %in% target_state_abbs
      owi_mask <- tolower(gb[[code_col_gb]]) %in% tolower(OWI_CODES)
      keep_mask_gb <- target_state_mask_gb & owi_mask
      
      cat(sprintf("  group_b: %s in 4 states, %s OWI, %s both\n",
                  format(sum(target_state_mask_gb), big.mark = ","),
                  format(sum(owi_mask), big.mark = ","),
                  format(sum(keep_mask_gb), big.mark = ",")))
      
      if (sum(keep_mask_gb) > 0) {
        collected[[length(collected) + 1]] <- data.table(
          ori = gb$ori[keep_mask_gb],
          offense_code = gb[[code_col_gb]][keep_mask_gb],
          year_inc = yr,
          source = "group_b"
        )
      }
    }
    rm(gb); gc(verbose = FALSE)
  } else {
    cat("  group_b: no file for", yr, "\n")
  }
}

# -----------------------------------------------------------------------------
# 5. Combine, classify, map to county
# -----------------------------------------------------------------------------
if (length(collected) == 0) stop("No data collected.")
ia <- rbindlist(collected, use.names = TRUE)
cat("\n=== Combined 4-state rows:", format(nrow(ia), big.mark = ","), "===\n")

ia[, offense_group := fcase(
  tolower(offense_code) %in% tolower(OWI_CODES), "owi",
  offense_code %in% DRUG_CODES,                   "drug",
  offense_code %in% PROPERTY_CODES,               "property",
  offense_code %in% VIOLENT_CODES,                "violent"
)]

ia <- merge(ia, leaic_4, by = "ori", all.x = TRUE)
n_unmatched <- sum(is.na(ia$county_fips))
cat("Unmatched ORIs:", n_unmatched,
    sprintf("(%.1f%%)\n", 100 * n_unmatched / nrow(ia)))

# -----------------------------------------------------------------------------
# 6. Aggregate
# -----------------------------------------------------------------------------
long <- ia[!is.na(county_fips) & !is.na(offense_group),
           .N, by = .(county_fips, state_name, year = year_inc, offense_group)]
wide <- dcast(long, county_fips + state_name + year ~ offense_group,
              value.var = "N", fill = 0)

for (col in c("owi", "drug", "property", "violent")) {
  if (!col %in% names(wide)) wide[, (col) := 0L]
}
wide <- wide[year %in% target_years]

dir.create(here::here("data", "interim"), showWarnings = FALSE, recursive = TRUE)
saveRDS(wide, here::here("data", "interim", "nibrs_county_year.rds"))

cat("\n=== DONE ===\n")
cat("Wrote nibrs_county_year.rds:", nrow(wide), "rows\n")
cat("  unique counties:", length(unique(wide$county_fips)), "\n")
cat("  state breakdown:\n")
print(wide[, .(n_county_years = .N,
               unique_counties = uniqueN(county_fips)),
           by = state_name])

cat("\n=== Outcome totals by state ===\n")
print(wide[, .(total_drug = sum(drug),
               total_owi  = sum(owi),
               total_prop = sum(property),
               total_viol = sum(violent)), by = state_name])
