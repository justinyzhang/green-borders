# ============================================================================
# R/04_cross_section_rebuild.R
# ----------------------------------------------------------------------------
# Rebuild Chapter 4 cross-section from raw IDFPR data.
# Verifies thesis1.Rmd's β = 2.895 result.
#
# Inputs:  /mnt/user-data/uploads/dispensaries_clean.csv (276 rows)
#          tidycensus ACS 2020 IL county population
#
# Outputs: out/ch4_cs_data.rds                 -- 46-county cross-section
#          out/ch4_cs_results.rds              -- 5 spec coefficients
#          out/ch4_cs_results.txt              -- human-readable table
#
# Sample: 2019-2024 cumulative new licenses (matches thesis1.Rmd panel cut)
# Specs:  (1) OLS bivariate    entry_per_100k ~ border
#         (2) OLS+Pop           entry_per_100k ~ border + log(population)
#         (3) OLS Trim          (1) excluding top 5 counties by total_entry
#         (4) PPML              total_entry ~ border + log(population)
#         (5) PPML Trim         (4) excluding top 5
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("dplyr", "readr", "fixest", "tidycensus")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(dplyr)
library(readr)
library(fixest)
library(tidycensus)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

# Helper: robust extraction (same as R_14/R_15 — no pvalue() dependency)
get_estimate <- function(m, var = "border") {
  if (inherits(m, "fixest")) {
    ct <- summary(m)$coeftable
  } else {
    ct <- summary(m)$coefficients
  }
  if (!var %in% rownames(ct)) return(list(beta = NA, se = NA, p = NA))
  # column names differ between lm and fixest; handle both
  p_col <- if ("Pr(>|t|)" %in% colnames(ct)) "Pr(>|t|)" else if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "p.value"
  list(
    beta = ct[var, "Estimate"],
    se   = ct[var, "Std. Error"],
    p    = ct[var, p_col]
  )
}

# ---- 1. Load dispensaries_clean.csv ---------------------------------------

disp_path <- "/mnt/user-data/uploads/dispensaries_clean.csv"
if (!file.exists(disp_path)) {
  # try local relative paths
  alt_paths <- c("dispensaries_clean.csv",
                 "data/raw/dispensaries_clean.csv",
                 "data/processed/dispensaries_clean.csv")
  for (p in alt_paths) {
    if (file.exists(p)) { disp_path <- p; break }
  }
}
stopifnot(file.exists(disp_path))

cat(sprintf("Reading: %s\n", disp_path))
disp <- read_csv(disp_path, show_col_types = FALSE)

cat(sprintf("Total dispensaries in file: %d\n", nrow(disp)))
cat(sprintf("Year range: %d to %d\n", min(disp$Year), max(disp$Year)))

# Normalize county names (handle DeKalb/Dekalb, Dupage/DuPage etc.)
disp <- disp %>%
  mutate(County = stringr::str_to_title(County)) %>%
  mutate(County = case_when(
    County == "Dekalb" ~ "DeKalb",
    County == "Dupage" ~ "DuPage",
    County == "Mchenry" ~ "McHenry",
    County == "Mclean" ~ "McLean",
    County == "Mcdonough" ~ "McDonough",
    County == "Lasalle" ~ "LaSalle",
    TRUE ~ County
  ))

# ---- 2. Restrict to 2019-2024 (matches thesis1.Rmd panel cut) -------------

disp_2024 <- disp %>% filter(Year >= 2019, Year <= 2024)
cat(sprintf("\n2019-2024 cumulative: %d dispensaries across %d counties\n",
            nrow(disp_2024), n_distinct(disp_2024$County)))

# Sanity check: Cook + Massac
cook_n <- sum(disp_2024$County == "Cook")
massac_n <- sum(disp_2024$County == "Massac")
cat(sprintf("Cook County: %d, Massac County: %d\n", cook_n, massac_n))

# ---- 3. Build county-level total_entry -----------------------------------

county_entry <- disp_2024 %>%
  count(County, name = "total_entry") %>%
  arrange(desc(total_entry))

cat("\n=== Top 10 counties by total_entry ===\n")
print(head(county_entry, 10))

# ---- 4. Border indicator --------------------------------------------------

# From dispensaries_clean.csv near_border column: take any() at county level
border_county <- disp_2024 %>%
  group_by(County) %>%
  summarise(
    border = as.integer(any(near_border)),
    avg_dist = mean(border_dist_mi, na.rm = TRUE)
  )

county_entry <- county_entry %>%
  left_join(border_county, by = "County")

cat(sprintf("\nBorder counties: %d / %d\n",
            sum(county_entry$border, na.rm = TRUE), nrow(county_entry)))

# ---- 5. Get ACS 2020 IL county population --------------------------------

cat("\nFetching ACS 2020 IL county population via tidycensus...\n")

# If user has Census API key set
tryCatch({
  if (Sys.getenv("CENSUS_API_KEY") == "") {
    cat("Note: CENSUS_API_KEY not set in environment.\n")
    cat("If this fails, get a free key from https://api.census.gov/data/key_signup.html\n")
    cat("Then run: census_api_key('YOUR_KEY', install = TRUE)\n")
    # try the key embedded in thesis1.Rmd
    api_key <- "06ccec4c76195625755962fe0462040884f4f1d1"
    cat("Trying embedded API key from thesis1.Rmd...\n")
    il_pop <- get_acs(
      geography = "county", state = "IL",
      variables = "B01003_001",
      year = 2020, survey = "acs5",
      key = api_key
    )
  } else {
    il_pop <- get_acs(
      geography = "county", state = "IL",
      variables = "B01003_001",
      year = 2020, survey = "acs5"
    )
  }
  
  il_pop <- il_pop %>%
    mutate(County = gsub(" County, Illinois", "", NAME)) %>%
    select(County, population = estimate)
  
  cat(sprintf("Loaded %d IL counties from ACS 2020\n", nrow(il_pop)))
  
}, error = function(e) {
  cat("\nERROR: tidycensus failed. Error message:\n", conditionMessage(e), "\n\n")
  cat("Fallback: hardcoded populations for key counties\n")
  # Hardcoded fallback for top counties
  il_pop <<- tibble::tribble(
    ~County, ~population,
    "Cook", 5275541,
    "DuPage", 932877,
    "Lake", 714342,
    "Will", 696355,
    "Kane", 516522,
    "McHenry", 310229,
    "Winnebago", 285350,
    "Madison", 264776,
    "St. Clair", 257400,
    "Champaign", 205865,
    "Sangamon", 196343,
    "Peoria", 179179,
    "Rock Island", 144477,
    "Tazewell", 130466,
    "Kankakee", 107502,
    "Macon", 103998,
    "McLean", 170628,
    "Vermilion", 74188,
    "Adams", 65737,
    "Massac", 14019,
    "Jackson", 56750,
    "Coles", 51188,
    "Effingham", 34242,
    "Jefferson", 37684,
    "Knox", 50321,
    "Logan", 28618,
    "Marion", 36899,
    "Montgomery", 28732,
    "Morgan", 33658,
    "Saline", 23768,
    "Union", 16653,
    "White", 13537,
    "Wabash", 11528,
    "Boone", 53606,
    "Jo Daviess", 21235,
    "Lawrence", 15678,
    "Alexander", 5761,
    "Fulton", 33609,
    "Grundy", 51992,
    "Kendall", 131869,
    "LaSalle", 109658,
    "Lee", 33828,
    "Livingston", 35648,
    "DeKalb", 100420,
    "Franklin", 38469
  )
})

# ---- 6. Merge to cross-section data --------------------------------------

cs_data <- county_entry %>%
  left_join(il_pop, by = "County") %>%
  filter(!is.na(population)) %>%
  mutate(entry_per_100k = total_entry / population * 100000)

cat(sprintf("\nCross-section built: %d counties\n", nrow(cs_data)))
cat(sprintf("Total entry summed: %d\n", sum(cs_data$total_entry)))
cat(sprintf("Counties with border==1: %d\n", sum(cs_data$border == 1)))
cat(sprintf("Counties with border==0: %d\n", sum(cs_data$border == 0)))

cat("\n=== Massac sanity check ===\n")
print(cs_data %>% filter(County == "Massac"))
cat("\n=== Cook sanity check ===\n")
print(cs_data %>% filter(County == "Cook"))

# ---- 7. Run 5 specifications --------------------------------------------

cat("\n=== Running 5 specifications ===\n\n")

# Top 5 counties to trim (for specs 3 and 5)
top5 <- cs_data %>% arrange(desc(total_entry)) %>% head(5) %>% pull(County)
cat(sprintf("Top 5 counties (trimmed in specs 3 & 5): %s\n\n",
            paste(top5, collapse = ", ")))

# Spec 1: OLS bivariate
m1 <- lm(entry_per_100k ~ border, data = cs_data)
e1 <- get_estimate(m1)
cat(sprintf("(1) OLS bivariate       β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e1$beta, e1$se, e1$p, nobs(m1)))

# Spec 2: OLS + log(pop)
m2 <- lm(entry_per_100k ~ border + log(population), data = cs_data)
e2 <- get_estimate(m2)
log_pop_e2 <- if ("log(population)" %in% rownames(summary(m2)$coefficients)) {
  summary(m2)$coefficients["log(population)", "Estimate"]
} else NA
cat(sprintf("(2) OLS+log(Pop)        β = %+.4f  SE = %.4f  p = %.4f  N = %d  [log(pop) coef = %.3f]\n",
            e2$beta, e2$se, e2$p, nobs(m2), log_pop_e2))

# Spec 3: OLS bivariate trim
cs_trim <- cs_data %>% filter(!County %in% top5)
m3 <- lm(entry_per_100k ~ border, data = cs_trim)
e3 <- get_estimate(m3)
cat(sprintf("(3) OLS Trim (drop top5) β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e3$beta, e3$se, e3$p, nobs(m3)))

# Spec 4: PPML
m4 <- fepois(total_entry ~ border + log(population), data = cs_data, vcov = "hetero")
e4 <- get_estimate(m4)
cat(sprintf("(4) PPML                β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e4$beta, e4$se, e4$p, nobs(m4)))

# Spec 5: PPML Trim
m5 <- fepois(total_entry ~ border + log(population), data = cs_trim, vcov = "hetero")
e5 <- get_estimate(m5)
cat(sprintf("(5) PPML Trim           β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e5$beta, e5$se, e5$p, nobs(m5)))

# ---- 8. Save outputs -----------------------------------------------------

results <- tibble::tibble(
  spec = c("(1) OLS bivariate", "(2) OLS+log(Pop)", "(3) OLS Trim",
           "(4) PPML", "(5) PPML Trim"),
  beta = c(e1$beta, e2$beta, e3$beta, e4$beta, e5$beta),
  se   = c(e1$se,   e2$se,   e3$se,   e4$se,   e5$se),
  p    = c(e1$p,    e2$p,    e3$p,    e4$p,    e5$p),
  n    = c(nobs(m1), nobs(m2), nobs(m3), nobs(m4), nobs(m5))
) %>%
  mutate(across(c(beta, se), ~round(., 4)),
         p = round(p, 4))

saveRDS(cs_data, "out/ch4_cs_data.rds")
saveRDS(results, "out/ch4_cs_results.rds")

sink("out/ch4_cs_results.txt")
cat("Chapter 4 cross-section verification\n")
cat("Source: dispensaries_clean.csv (2019-2024 cumulative)\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Total dispensaries: %d\n", sum(cs_data$total_entry)))
cat(sprintf("Counties with positive entry: %d\n", nrow(cs_data)))
cat(sprintf("Border counties: %d\n", sum(cs_data$border == 1)))
cat(sprintf("Cook: %d, Massac: %d\n",
            cs_data$total_entry[cs_data$County == "Cook"],
            cs_data$total_entry[cs_data$County == "Massac"]))
cat(strrep("-", 70), "\n", sep = "")
print(results)
sink()

cat("\n=== Saved ===\n")
cat("  out/ch4_cs_data.rds\n")
cat("  out/ch4_cs_results.rds\n")
cat("  out/ch4_cs_results.txt\n")

# ---- 9. Verification vs Defense Speech v2 -------------------------------

cat("\n", strrep("=", 70), "\n", sep = "")
cat("VERIFICATION vs Defense Speech v2 (April 2026)\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Speech claim: 113 dispensaries, 46 counties, Cook ?, Massac 3\n"))
cat(sprintf("Reality:      %d dispensaries, %d counties, Cook %d, Massac %d\n",
            sum(cs_data$total_entry), nrow(cs_data),
            cs_data$total_entry[cs_data$County == "Cook"],
            cs_data$total_entry[cs_data$County == "Massac"]))
cat(sprintf("\nSpeech claim: β(OLS+Pop) = 2.895, p = 0.003\n"))
cat(sprintf("Reality:      β(OLS+Pop) = %.3f, p = %.3f\n", e2$beta, e2$p))

if (abs(e2$beta - 2.895) < 0.01) {
  cat("✓ β MATCHES defense speech (within rounding)\n")
} else {
  cat(sprintf("⚠️ β DIFFERS by %.4f — possible reasons:\n", abs(e2$beta - 2.895)))
  cat("   - Different border indicator definition\n")
  cat("   - Different sample (population filter, missing data)\n")
  cat("   - Different ACS year\n")
}

cat("\nDone. Paste this output to your Claude session.\n")
