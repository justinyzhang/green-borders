# =============================================================================
# 02_acs.R   (4-state version: IA + WI + IN + MO)
#
# Pull ACS 5-year estimates 2015-2022 for all counties in 4 prohibition states
# via direct Census API (no tidycensus dependency).
#
# Output:
#   data/interim/acs_iowa.rds   (kept name for downstream compatibility)
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(data.table)
  library(here)
})

CENSUS_API_KEY <- Sys.getenv("CENSUS_API_KEY", unset = "")
if (!nchar(CENSUS_API_KEY)) {
  stop("CENSUS_API_KEY environment variable not set. ",
       "Get a free key at https://api.census.gov/data/key_signup.html ",
       "and add it to your .Renviron file.")
}

acs_vars <- c(
  pop_total = "B01003_001E",
  med_inc   = "B19013_001E",
  pov_total = "B17001_002E",
  pop_white = "B02001_002E",
  pop_hisp  = "B03002_012E"
)

target_years <- 2015:2022
target_state_fips <- c("19", "55", "18", "29")   # IA, WI, IN, MO

fetch_acs_year_state <- function(yr, state_fips) {
  vars <- paste(c("NAME", unname(acs_vars)), collapse = ",")
  url <- sprintf(
    "https://api.census.gov/data/%d/acs/acs5?get=%s&for=county:*&in=state:%s&key=%s",
    yr, vars, state_fips, CENSUS_API_KEY
  )
  resp <- httr::GET(url)
  if (httr::status_code(resp) != 200) {
    stop(sprintf("Census API status %d for year %d state %s",
                 httr::status_code(resp), yr, state_fips))
  }
  raw <- httr::content(resp, "text", encoding = "UTF-8")
  parsed <- jsonlite::fromJSON(raw)
  hdr <- parsed[1, ]; body <- parsed[-1, , drop = FALSE]
  df <- as.data.frame(body, stringsAsFactors = FALSE)
  names(df) <- hdr
  for (v in unname(acs_vars)) {
    if (v %in% names(df)) df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
  }
  df$GEOID <- paste0(df$state, df$county)
  df$year <- yr
  as.data.table(df)
}

cat("Fetching ACS 5-year for 4 states x 8 years = 32 API calls ...\n")
all <- list()
for (yr in target_years) {
  for (st in target_state_fips) {
    cat(sprintf("  %d state %s ...\n", yr, st))
    all[[length(all) + 1]] <- fetch_acs_year_state(yr, st)
  }
}
acs <- rbindlist(all, use.names = TRUE, fill = TRUE)

for (i in seq_along(acs_vars)) {
  api_name <- unname(acs_vars[i]); friendly <- names(acs_vars)[i]
  if (api_name %in% names(acs)) setnames(acs, api_name, friendly)
}

acs_slim <- acs[, .(
  county_fips = GEOID,
  year,
  pop_total, med_inc, pov_total, pop_white, pop_hisp
)]

dir.create(here::here("data", "interim"), showWarnings = FALSE, recursive = TRUE)
saveRDS(acs_slim, here::here("data", "interim", "acs_iowa.rds"))

cat("\n=== DONE ===\n")
cat("Wrote acs_iowa.rds:", nrow(acs_slim), "rows\n")
cat("  counties:", length(unique(acs_slim$county_fips)), "\n")
cat("  years:   ", paste(range(acs_slim$year), collapse = " - "), "\n")
