# Claude Code Task Examples

Reference patterns for issuing tasks to Claude Code. Each task should:

1. Reference a **cluster number** from `docs/research_design.md`.
2. Name the script files to read or modify.
3. Specify inputs (file paths) and outputs (file paths or stdout format).
4. State acceptance criteria (test that should pass, table column count, etc.).

## Tier 0 — Skeleton tasks (do these first, before any analysis)

### T0.1 Build project skeleton

> Build the project skeleton. Create the directory structure described in `CLAUDE.md`. Initialize an `renv` environment locked to the package set: `data.table`, `here`, `sandwich`, `lmtest`, `ggplot2`, `httr`, `jsonlite`, `pdftools`, `stringr`, `sf`, `tigris`, `fixest`, `did`, `bacondecomp`. Do **not** add `arrow` or `tidycensus` — those have rlang incompatibilities (see `CLAUDE.md`). Write a placeholder `R/00_master.R` that sources scripts `00_build_metadata.R` through `09_sumstats.R` in order, with `if (file.exists(...))` guards so missing scripts don't crash the master.

### T0.2 Add gitignore

> Create `.gitignore`. Ignore: `data/raw/` (everything), `data/interim/` (everything), `*.Rhistory`, `.Rproj.user/`, `*.Rdata`, `renv/library/`, `renv/python/`, `renv/staging/`, `.Renviron`, `*.aux`, `*.log`, `*.bbl`, `*.blg`, `*.toc`. Do **not** ignore `data/processed/` (the small final panel is fine to track) or `output/`.

## Tier 1 — Utility functions

### T1.1 Spatial aggregation helper (cluster 11)

> Read `docs/research_design.md` cluster 11. Implement a function `compute_radius_population(county_fips, radius_miles, state_subset)` in `R/utils/spatial_aggregation.R` that, given an origin county FIPS, returns the total population from counties whose centroids lie within `radius_miles` of the origin centroid, optionally restricted to a vector of 2-digit state FIPS codes (`state_subset`). Use `data/raw/county_centroids.csv` (columns: `county_fips`, `lat`, `lon`) and `data/raw/county_pop_2020.csv` (columns: `county_fips`, `pop_2020`). Use the haversine formula; don't depend on `sf` for this function. Add unit tests in `tests/test_spatial_aggregation.R` covering: (a) origin within 15 miles of itself returns own population, (b) Cook County IL with `state_subset = c("17")` returns within an order of magnitude of 5.2M, (c) a Lake Michigan-coastal county returns a stable result despite missing centroids in the lake.

### T1.2 Iowa DOT PDF parser

> Implement `parse_iadot_pdf(pdf_path, year_range, outcome_label, iowa_counties)` in `R/utils/pdf_parser.R`. Input is a single Iowa DOT statistics PDF (county × year wide table). Use `pdftools::pdf_text()` to extract text, then regex-match each line against the 99-county name vector (all-uppercase, multi-word names like "BLACK HAWK", "CERRO GORDO", "DES MOINES", "PALO ALTO", "VAN BUREN"). Return a `data.table` with columns `name_upper`, `year`, `value`, `outcome`. If a CSV file exists at `sub("\\.pdf$", ".csv", pdf_path)`, read it directly and skip PDF parsing (manual fallback). Unit-test against a small fixture file with 3 counties × 5 years.

### T1.3 TWFE without fixest

> Create `R/utils/twfe_no_fixest.R` exporting two functions: `fit_twfe(formula_rhs, data, cluster_col, weights = NULL)` and `fit_event_study(y_col, event_var, ref_k, ...)`. Both use base `lm()` plus `sandwich::vcovCL(type = "HC1")` plus `lmtest::coeftest`. Return a list with `beta`, `se`, `p`, `ci_lo`, `ci_hi`, `n_obs`, `n_cluster`. Match the API of `fixest::feols` for the parts the rest of the pipeline uses (so a future swap is mechanical). Tests in `tests/test_twfe.R` should compare against `fixest::feols` output on a synthetic 50-county × 8-year panel when fixest is installed (use `testthat::skip_if_not_installed("fixest")`).

## Tier 2 — Pipeline scripts (one at a time)

### T2.1 NIBRS loader

> Read `docs/research_design.md` cluster 4 and the comments at the top of `R/01_nibrs_load.R`. Verify the script: (a) handles all 4 prohibition states (IA, IN, WI, MO) via the `target_state_names` vector, (b) correctly auto-detects the offense code column in the Group B segment (it may be `ucr_arrest_offense_code` or `arrest_offense_code`), (c) writes `data/interim/nibrs_county_year.rds` with one row per county-year and columns `county_fips`, `state_name`, `year`, `drug`, `owi`, `property`, `violent`. Print a summary table showing rows per state. Do not modify the offense-code synonym lists.

### T2.2 ACS pull

> Verify `R/02_acs.R` pulls ACS 5-year estimates for 4 states (FIPS 18, 19, 29, 55) for years 2015-2022. Use direct HTTP to the Census API (not `tidycensus`, which has rlang issues). Read the API key from `Sys.getenv("CENSUS_API_KEY")`; print a warning if unset. Output `data/interim/acs_iowa.rds` (filename kept for downstream compatibility) with one row per county-year and columns `county_fips`, `year`, `pop_total`, `med_inc`, `pov_total`, `pop_white`, `pop_hisp`. Total expected rows: ~378 counties × 8 years = ~3024.

### T2.3 Panel assembler

> Read `docs/research_design.md` clusters 1 + 2 + 4. Verify `R/03_panel.R` builds a balanced panel by merging `data/processed/iowa_county_metadata.csv` (county metadata, treated assignment), `data/interim/nibrs_county_year.rds` (NIBRS counts), and `data/interim/acs_iowa.rds` (demographics). The output `data/processed/panel_county_year.rds` should have exactly 263 × 8 = 2104 rows for the main DiD sample (treatment %in% c("IL_border", "interior")), plus an additional ~115 × 8 rows for Missouri held-out placebo. Verify `treat = is_il_border * post` is correctly computed and equals 1 for exactly 81 obs (27 treated counties × 3 post years).

## Tier 3 — Estimation

### T3.1 Baseline DiD (cluster 5)

> Read `docs/research_design.md` cluster 5. In `R/05_twfe.R`, fit two specifications for each of four outcomes (`ln_drug`, `ln_owi`, `ln_property`, `ln_violent`):
> - Spec A: county FE + state-by-year FE.
> - Spec B: county FE + year FE.
> Sample: `treatment %in% c("IL_border", "interior")` (excluding Missouri).
> Cluster SE at county level via `sandwich::vcovCL(type = "HC1")`.
> Write `output/tables/tab1_baseline_coefs.csv` with columns: `outcome`, `spec`, `beta`, `se`, `p_value`, `ci_lo`, `ci_hi`, `n_obs`.
> Write `output/tables/tab1_baseline.tex` as a 4-column booktabs LaTeX table reporting Spec A only.

### T3.2 Event study (cluster 7)

> Read `docs/research_design.md` cluster 7. In `R/06_eventstudy.R`, fit event-study specifications for all 4 outcomes. Window k ∈ {-5..2}, reference k = -1. Use the `R/utils/twfe_no_fixest.R` event-study helper. For each outcome, write a PDF figure with `ggplot2`: errorbars + points, vertical reference line at k = -0.5, horizontal at 0. Write all coefficients to `output/tables/fig3_eventstudy_coefs.csv` (long format: `outcome`, `k`, `est`, `se`, `ci_lo`, `ci_hi`).

### T3.3 Robustness (cluster 8)

> Read `docs/research_design.md` cluster 8. In `R/07_robust.R`, run the 7 robustness specs on `ln_drug` only. Output `output/tables/tab2_robust_coefs.csv` (one row per spec) and `output/tables/tab2_robust.tex` (5 main specs as columns; per-state results go in a separate table).

### T3.4 Gravity entry regression (cluster 6)

> Read `docs/research_design.md` cluster 6. Implement `R/04_gravity_regression.R` from scratch. Inputs: `data/processed/il_county_entry.csv` (operational dispensary count per IL county) and the spatial aggregation helper from T1.1 to compute `Pop_IL_X` and `Pop_PROH_X` at X = 15, 50, 100 miles. Estimator: `glm(family = quasipoisson(link = "log"))` if `fixest::fepois` fails to load. Report three columns in `output/tables/gravity_main.tex`:
> 1. Binary `Border_15` dummy only.
> 2. `Border_15 * log(Pop_PROH_15)` interaction.
> 3. `log(Pop_PROH_15)` in levels.
> Cluster SE at county. Include each county's own `log(Pop_IL_15)` as a control in all columns.

### T3.5 Iowa DOT registry DiD (cluster 4b)

> Read `docs/research_design.md` cluster 4b. Verify `R/08_iowadot_did.R` runs the same Iowa-only DiD spec from cluster 5 on the three Iowa DOT outcomes: OWI revocations, speeding convictions, traffic fatalities. Sample: 99 Iowa counties × 8 years. Output `output/tables/tab3_iowadot_baseline.csv` with one row per outcome. Output `output/figures/fig4_iowadot_eventstudy.pdf` with a 3-panel facet (one panel per outcome).

## Tier 4 — Reporting

### T4.1 Summary stats table

> Implement `R/09_sumstats.R`. Read `data/processed/panel_county_year.rds`. Generate a balance table: for each variable in {`pop_total`, `med_inc`, `pov_rate`, `drug`, `owi`, `property`, `violent`}, report mean and SD in the pre-period (2015-2019), separately for treated (border = 1) and control (border = 0), plus a t-test of the difference. Write `output/tables/tab0_sumstats.tex`.

### T4.2 Master orchestration

> Update `R/00_master.R` to source `00_build_metadata.R` through `09_sumstats.R` in order. Each step should be guarded by `tryCatch` so a failure in one step doesn't halt the others. Print elapsed wall time per step. Print a final summary table of file sizes for everything in `output/`.

## Anti-patterns to avoid

- Do **not** generate Claude Code tasks that say "make the regression significant" or "search for an outcome that gives p < 0.05". That's specification-mining. Tasks should reference clusters, not p-values.
- Do **not** ask Claude Code to invent specifications. The clusters are fixed; new specifications go into a new numbered cluster after user approval.
- Do **not** combine multiple unrelated changes in one task. One file, one cluster, one acceptance criterion at a time.
