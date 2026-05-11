# =============================================================================
# 01c_iowadot_load.R   (Iowa DOT OWI revocations + speeding citations panel)
#
# Pull 3 Iowa DOT county-year statistics:
#   1. OWI revocations by county (2005-2024)
#   2. Speeding convictions by county (2013-2022)
#   3. Yearly traffic fatalities by county (2015-2024)
#
# All three are clean state-government registry data, not subject to NIBRS
# reporting coverage rotation issues that plagued the NIBRS analysis.
#
# Strategy:
#   - Download PDFs from iowadot.gov/mvd/FactsandStats
#   - Parse with pdftools::pdf_text + regex
#   - Manual CSV fallback if parsing fails
#
# Outputs:
#   data/raw/iowadot/owi_revocations.pdf
#   data/raw/iowadot/speeding_convictions.pdf
#   data/raw/iowadot/yearly_fatalities.pdf
#   data/interim/iowadot_county_year.rds   (long panel)
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(pdftools)
  library(data.table)
  library(stringr)
  library(here)
})

dir.create(here::here("data", "raw", "iowadot"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("data", "interim"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. URLs (verified from iowadot.gov/mvd/FactsandStats)
# -----------------------------------------------------------------------------
# These media IDs are the working URLs as of search verification.
# If a URL 404s, browse iowadot.gov/mvd/FactsandStats and right-click the
# right link to get the current /media/NNNN/download URL.
urls <- list(
  owi      = list(url = "https://iowadot.gov/media/7249/download",
                  file = "owi_revocations.pdf",
                  outcome = "owi_revocation",
                  year_range = 2005:2024),
  speeding = list(url = "https://iowadot.gov/mvd/FactsandStats",  # fallback: find on page
                  file = "speeding_convictions.pdf",
                  outcome = "speeding_conviction",
                  year_range = 2013:2022),
  fatality = list(url = "https://iowadot.gov/mvd/FactsandStats",  # fallback: find on page
                  file = "yearly_fatalities.pdf",
                  outcome = "fatality",
                  year_range = 2015:2024)
)

# -----------------------------------------------------------------------------
# 2. Download
# -----------------------------------------------------------------------------
for (nm in names(urls)) {
  ent <- urls[[nm]]
  out_path <- here::here("data", "raw", "iowadot", ent$file)
  if (file.exists(out_path) && file.size(out_path) > 1000) {
    cat("Skip:", ent$file, "already downloaded\n")
    next
  }
  cat("Downloading", nm, "from", ent$url, "...\n")
  resp <- try(httr::GET(ent$url,
                        httr::write_disk(out_path, overwrite = TRUE),
                        httr::user_agent("Mozilla/5.0 (research)")),
              silent = TRUE)
  if (inherits(resp, "try-error") || httr::status_code(resp) != 200) {
    cat("  FAIL: download failed for", nm, "\n")
    cat("  ACTION: manually save", ent$url, "as", out_path, "\n")
    cat("  OR: convert PDF to CSV manually, place as",
        sub("\\.pdf$", ".csv", out_path), "\n")
  } else {
    cat("  OK:", file.size(out_path), "bytes\n")
  }
}

# -----------------------------------------------------------------------------
# 3. Iowa 99 county names (alphabetical, official spelling matching DOT)
# -----------------------------------------------------------------------------
iowa_counties <- c(
  "ADAIR", "ADAMS", "ALLAMAKEE", "APPANOOSE", "AUDUBON", "BENTON",
  "BLACK HAWK", "BOONE", "BREMER", "BUCHANAN", "BUENA VISTA", "BUTLER",
  "CALHOUN", "CARROLL", "CASS", "CEDAR", "CERRO GORDO", "CHEROKEE",
  "CHICKASAW", "CLARKE", "CLAY", "CLAYTON", "CLINTON", "CRAWFORD",
  "DALLAS", "DAVIS", "DECATUR", "DELAWARE", "DES MOINES", "DICKINSON",
  "DUBUQUE", "EMMET", "FAYETTE", "FLOYD", "FRANKLIN", "FREMONT",
  "GREENE", "GRUNDY", "GUTHRIE", "HAMILTON", "HANCOCK", "HARDIN",
  "HARRISON", "HENRY", "HOWARD", "HUMBOLDT", "IDA", "IOWA",
  "JACKSON", "JASPER", "JEFFERSON", "JOHNSON", "JONES", "KEOKUK",
  "KOSSUTH", "LEE", "LINN", "LOUISA", "LUCAS", "LYON", "MADISON",
  "MAHASKA", "MARION", "MARSHALL", "MILLS", "MITCHELL", "MONONA",
  "MONROE", "MONTGOMERY", "MUSCATINE", "OBRIEN", "OSCEOLA", "PAGE",
  "PALO ALTO", "PLYMOUTH", "POCAHONTAS", "POLK", "POTTAWATTAMIE",
  "POWESHIEK", "RINGGOLD", "SAC", "SCOTT", "SHELBY", "SIOUX", "STORY",
  "TAMA", "TAYLOR", "UNION", "VAN BUREN", "WAPELLO", "WARREN",
  "WASHINGTON", "WAYNE", "WEBSTER", "WINNEBAGO", "WINNESHIEK",
  "WOODBURY", "WORTH", "WRIGHT"
)
stopifnot(length(iowa_counties) == 99)

# Map county names to 5-digit FIPS codes (state 19)
# Order matches alphabetical FIPS order, which is exactly the order above
iowa_fips <- sprintf("19%03d", seq(1, 197, by = 2))
county_fips_map <- data.table(name_upper = iowa_counties, county_fips = iowa_fips)

# -----------------------------------------------------------------------------
# 4. PDF parser - generic wide table where rows = county, cols = years
# -----------------------------------------------------------------------------
parse_iadot_pdf <- function(pdf_path, year_range, outcome_label) {
  if (!file.exists(pdf_path)) {
    csv_path <- sub("\\.pdf$", ".csv", pdf_path)
    if (file.exists(csv_path)) {
      cat("  Using manual CSV at", csv_path, "\n")
      df <- fread(csv_path)
      return(df)
    }
    cat("  PDF not found and no manual CSV. Skipping.\n")
    return(NULL)
  }
  
  txt <- tryCatch(pdftools::pdf_text(pdf_path),
                  error = function(e) {cat("  ERROR reading PDF:", conditionMessage(e), "\n"); NULL})
  if (is.null(txt)) return(NULL)
  
  rows <- list()
  for (page_text in txt) {
    lines <- strsplit(page_text, "\n")[[1]]
    for (ln in lines) {
      ln <- str_trim(ln)
      if (nchar(ln) == 0) next
      # Try to match: <county name> <number> <number> ... <number>
      # County name can be 1-2 words. We try longest match first.
      for (cty in iowa_counties) {
        # Allow case-insensitive match with optional whitespace
        pattern <- paste0("^", gsub(" ", "\\\\s+", cty), "\\b")
        if (str_detect(toupper(ln), pattern)) {
          # Strip the county name from the line
          rest <- str_trim(sub(pattern, "", toupper(ln), perl = TRUE))
          # Extract all numbers (handle commas in thousands)
          nums <- as.integer(gsub(",", "",
                                   str_extract_all(rest, "[0-9,]+")[[1]]))
          nums <- nums[!is.na(nums)]
          if (length(nums) >= length(year_range)) {
            # Take exactly length(year_range) - assume left-aligned
            nums <- nums[seq_along(year_range)]
            rows[[length(rows) + 1]] <- data.table(
              name_upper = cty,
              year = year_range,
              value = nums
            )
          }
          break  # found county for this line
        }
      }
    }
  }
  
  if (length(rows) == 0) {
    cat("  WARN: parser found 0 county rows in", pdf_path, "\n")
    cat("  ACTION: open PDF in Excel/Acrobat, save as CSV at",
        sub("\\.pdf$", ".csv", pdf_path), "\n")
    return(NULL)
  }
  
  dt <- unique(rbindlist(rows), by = c("name_upper", "year"))
  dt[, outcome := outcome_label]
  cat(sprintf("  Parsed: %d counties x %d years from %s\n",
              uniqueN(dt$name_upper), uniqueN(dt$year), basename(pdf_path)))
  dt
}

# -----------------------------------------------------------------------------
# 5. Parse all 3 PDFs
# -----------------------------------------------------------------------------
all_data <- list()
for (nm in names(urls)) {
  ent <- urls[[nm]]
  cat("\nParsing", nm, "...\n")
  parsed <- parse_iadot_pdf(
    here::here("data", "raw", "iowadot", ent$file),
    ent$year_range, ent$outcome
  )
  if (!is.null(parsed)) all_data[[nm]] <- parsed
}

if (length(all_data) == 0) {
  stop("No data parsed. Manual CSV conversion required - see notes above.")
}

# -----------------------------------------------------------------------------
# 6. Combine, attach county FIPS, reshape to wide county-year panel
# -----------------------------------------------------------------------------
long <- rbindlist(all_data, use.names = TRUE, fill = TRUE)
long <- merge(long, county_fips_map, by = "name_upper", all.x = TRUE)

n_unmatched <- sum(is.na(long$county_fips))
if (n_unmatched > 0) {
  cat("WARN:", n_unmatched, "rows have unmatched county names\n")
  print(unique(long[is.na(county_fips), .(name_upper)]))
}
long <- long[!is.na(county_fips)]

wide <- dcast(long, county_fips + year ~ outcome, value.var = "value", fill = NA)

# Restrict to our analysis window 2015-2022 (same as NIBRS)
wide <- wide[year %in% 2015:2022]

saveRDS(wide, here::here("data", "interim", "iowadot_county_year.rds"))

cat("\n=== DONE ===\n")
cat("Wrote iowadot_county_year.rds:", nrow(wide), "rows\n")
cat("  unique counties:", length(unique(wide$county_fips)), "\n")
cat("  outcomes present:", paste(setdiff(names(wide), c("county_fips","year")), collapse = ", "), "\n")

cat("\n=== Scott County (19163, IL-border, Quad Cities) ===\n")
print(wide[county_fips == "19163"][order(year)])

cat("\n=== Polk County (19153, interior, Des Moines) ===\n")
print(wide[county_fips == "19153"][order(year)])

cat("\n=== Outcome totals across 99 counties, 2015-2022 ===\n")
for (col in setdiff(names(wide), c("county_fips","year"))) {
  cat(sprintf("  %s: total = %s, mean per county-year = %.1f\n",
              col, format(sum(wide[[col]], na.rm = TRUE), big.mark = ","),
              mean(wide[[col]], na.rm = TRUE)))
}
