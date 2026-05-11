# =============================================================================
# 00_master.R
#
# Master script. Sources analysis scripts in order, with tryCatch guards so a
# failure in one step does not halt downstream steps that have other inputs.
#
# Run with: source("R/00_master.R")
# =============================================================================

suppressPackageStartupMessages({
  library(here)
})

scripts <- c(
  "R/00_build_metadata.R",
  "R/01_nibrs_load.R",
  "R/01c_iowadot_load.R",
  "R/02_acs.R",
  "R/03_panel.R",
  "R/04_raw_plot.R",
  "R/04_gravity_regression.R",
  "R/05_twfe.R",
  "R/06_eventstudy.R",
  "R/07_robust.R",
  "R/08_iowadot_did.R",
  "R/09_sumstats.R"
)

run_one <- function(path) {
  if (!file.exists(here::here(path))) {
    cat(sprintf("[SKIP] %s (not present)\n", path))
    return(invisible(NULL))
  }
  cat(sprintf("[RUN]  %s\n", path))
  t0 <- Sys.time()
  out <- tryCatch(
    source(here::here(path), local = new.env()),
    error = function(e) {
      cat(sprintf("[FAIL] %s -> %s\n", path, conditionMessage(e)))
      NULL
    }
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("       %.1fs\n", dt))
}

for (s in scripts) run_one(s)

cat("\n=== Output summary ===\n")
out_dir <- here::here("output")
if (dir.exists(out_dir)) {
  files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  if (length(files) > 0) {
    info <- file.info(files)
    info$path <- files
    info <- info[order(-info$size), c("path", "size")]
    info$kb <- round(info$size / 1024, 1)
    print(head(info[, c("path", "kb")], 30))
  } else {
    cat("(empty)\n")
  }
}
