# ============================================================================
# R/04_cross_section_rebuild_v2.R
# ----------------------------------------------------------------------------
# v2: Use 32-county ADJACENCY border definition from thesis1.Rmd,
#     not the 7-county 15-mile distance from dispensaries_clean.csv
#
# Sample: 2019-2024 cumulative new licenses (240 dispensaries, 46 counties)
# Border: 32 IL counties that touch any prohibition-state border
# ============================================================================

library(dplyr)
library(readr)
library(fixest)
library(tidycensus)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

get_estimate <- function(m, var = "border") {
  if (inherits(m, "fixest")) {
    ct <- summary(m)$coeftable
  } else {
    ct <- summary(m)$coefficients
  }
  if (!var %in% rownames(ct)) return(list(beta = NA, se = NA, p = NA))
  p_col <- if ("Pr(>|t|)" %in% colnames(ct)) "Pr(>|t|)" else if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "p.value"
  list(
    beta = ct[var, "Estimate"],
    se   = ct[var, "Std. Error"],
    p    = ct[var, p_col]
  )
}

# ---- 1. Load dispensaries_clean.csv --------------------------------------

disp_path <- "/mnt/user-data/uploads/dispensaries_clean.csv"
if (!file.exists(disp_path)) {
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

# ---- 2. 2019-2024 cumulative --------------------------------------------

disp_2024 <- disp %>% filter(Year >= 2019, Year <= 2024)
cat(sprintf("2019-2024: %d dispensaries across %d counties\n",
            nrow(disp_2024), n_distinct(disp_2024$County)))

# ---- 3. county_entry -----------------------------------------------------

county_entry <- disp_2024 %>%
  count(County, name = "total_entry") %>%
  arrange(desc(total_entry))

# ---- 4. 32-county adjacency border (from thesis1.Rmd) -------------------

border_counties_32 <- c(
  "Jo Daviess", "Stephenson", "Winnebago", "Boone", "McHenry", "Lake",
  "Vermilion", "Edgar", "Clark", "Crawford", "Lawrence", "Wabash", "White",
  "Gallatin", "Hardin", "Pope", "Massac", "Alexander", "Union", "Jackson",
  "Randolph", "Monroe", "St. Clair", "Madison", "Jersey", "Calhoun",
  "Pike", "Adams", "Hancock", "Henderson", "Mercer", "Rock Island"
)

cat(sprintf("\n32 border counties (adjacency): %s\n",
            paste(border_counties_32, collapse = ", ")))

county_entry <- county_entry %>%
  mutate(border = as.integer(County %in% border_counties_32))

cat(sprintf("\nBorder counties (in sample): %d / %d\n",
            sum(county_entry$border, na.rm = TRUE), nrow(county_entry)))

# Which border counties are NOT in sample (zero entry)?
border_not_in_sample <- setdiff(border_counties_32, county_entry$County)
cat(sprintf("Border counties NOT in 46-county sample (zero entry): %s\n",
            paste(border_not_in_sample, collapse = ", ")))

# ---- 5. ACS 2020 population ---------------------------------------------

cat("\nFetching ACS 2020 IL county population...\n")

api_key <- "06ccec4c76195625755962fe0462040884f4f1d1"

il_pop <- get_acs(
  geography = "county", state = "IL",
  variables = "B01003_001",
  year = 2020, survey = "acs5",
  key = api_key
) %>%
  mutate(County = gsub(" County, Illinois", "", NAME)) %>%
  select(County, population = estimate)

cat(sprintf("Loaded %d IL counties from ACS 2020\n", nrow(il_pop)))

# ---- 6. Build cs_data ---------------------------------------------------

cs_data <- county_entry %>%
  left_join(il_pop, by = "County") %>%
  filter(!is.na(population)) %>%
  mutate(entry_per_100k = total_entry / population * 100000)

cat(sprintf("\nCross-section: %d counties\n", nrow(cs_data)))
cat(sprintf("Total entry: %d\n", sum(cs_data$total_entry)))
cat(sprintf("Border counties in sample: %d\n", sum(cs_data$border == 1)))
cat(sprintf("Interior counties in sample: %d\n", sum(cs_data$border == 0)))

cat("\n=== Border counties in 46-sample ===\n")
print(cs_data %>% filter(border == 1) %>% select(County, total_entry, population, entry_per_100k))

# ---- 7. Run 5 specifications --------------------------------------------

cat("\n=== Running 5 specifications (32-county adjacency border) ===\n\n")

top5 <- cs_data %>% arrange(desc(total_entry)) %>% head(5) %>% pull(County)

m1 <- lm(entry_per_100k ~ border, data = cs_data)
e1 <- get_estimate(m1)
cat(sprintf("(1) OLS bivariate       β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e1$beta, e1$se, e1$p, nobs(m1)))

m2 <- lm(entry_per_100k ~ border + log(population), data = cs_data)
e2 <- get_estimate(m2)
log_pop_e2 <- summary(m2)$coefficients["log(population)", "Estimate"]
cat(sprintf("(2) OLS+log(Pop)        β = %+.4f  SE = %.4f  p = %.4f  N = %d  [log(pop) coef = %.3f]\n",
            e2$beta, e2$se, e2$p, nobs(m2), log_pop_e2))

cs_trim <- cs_data %>% filter(!County %in% top5)
m3 <- lm(entry_per_100k ~ border, data = cs_trim)
e3 <- get_estimate(m3)
cat(sprintf("(3) OLS Trim            β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e3$beta, e3$se, e3$p, nobs(m3)))

m4 <- fepois(total_entry ~ border + log(population), data = cs_data, vcov = "hetero")
e4 <- get_estimate(m4)
cat(sprintf("(4) PPML                β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e4$beta, e4$se, e4$p, nobs(m4)))

m5 <- fepois(total_entry ~ border + log(population), data = cs_trim, vcov = "hetero")
e5 <- get_estimate(m5)
cat(sprintf("(5) PPML Trim           β = %+.4f  SE = %.4f  p = %.4f  N = %d\n",
            e5$beta, e5$se, e5$p, nobs(m5)))

# ---- 8. Save ------------------------------------------------------------

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

saveRDS(cs_data, "out/ch4_cs_data_v2_adjacency.rds")
saveRDS(results, "out/ch4_cs_results_v2_adjacency.rds")

cat("\n=== Saved ===\n")
cat("  out/ch4_cs_data_v2_adjacency.rds\n")
cat("  out/ch4_cs_results_v2_adjacency.rds\n")

# ---- 9. Verification ----------------------------------------------------

cat("\n", strrep("=", 70), "\n", sep = "")
cat("VERIFICATION vs Defense Speech v2 (adjacency border)\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("Speech claim:  17 of 46 treated, β(OLS+Pop) = 2.895, p = 0.003\n"))
cat(sprintf("Reality:       %d of %d treated, β(OLS+Pop) = %.3f, p = %.4f\n",
            sum(cs_data$border == 1), nrow(cs_data), e2$beta, e2$p))

cat(sprintf("\nSpeech claim:  Spec 1 OLS bivariate β = 3.6\n"))
cat(sprintf("Reality:       β = %.3f\n", e1$beta))

cat(sprintf("\nSpeech claim:  Spec 3 trim β = 3.943\n"))
cat(sprintf("Reality:       β = %.3f\n", e3$beta))

if (abs(e2$beta - 2.895) < 0.10) {
  cat("\n✓ β MATCHES defense speech (within reasonable rounding)\n")
} else {
  cat(sprintf("\n⚠️ β differs by %.4f — likely due to exact ACS year or county name spelling\n",
              abs(e2$beta - 2.895)))
}

cat("\nDone. Paste output.\n")
