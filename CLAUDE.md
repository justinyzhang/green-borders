# CLAUDE.md

This file gives Claude Code the project-level context for the `green-borders` repository. Read it before any task.

## Project

Illinois recreational cannabis legalization (Jan 2020), cross-border dispensary entry, and downstream crime spillovers in IL-bordering counties of prohibition states.

- **Primary outcome** (job-market paper): NIBRS crime DiD in IL-border counties across Iowa, Indiana, and Wisconsin (with Missouri held out as a 2023 placebo).
- **Supporting outcome** (Chapter 4): PPML gravity entry regression on 46 Illinois counties with cross-border population pulls.
- **Defense date**: May 20, 2026.

## Conventions

- **Language**: R 4.5.x. Stata files only where explicitly noted.
- **Estimation**: `fixest` for fixed-effects regressions whenever the `rlang` version allows. If `fixest` fails to load (rlang < 1.1.7), fall back to base `lm()` + `sandwich::vcovCL` for cluster-robust SEs. A working fallback already exists in `R/utils/twfe_no_fixest.R`; reuse it.
- **Spatial**: `sf` for geometry, `tigris` for county polygons when available. Centroids and distances always in miles (haversine).
- **Data manipulation**: `data.table` is preferred for the panel construction scripts (`01_*` through `03_*`) because the NIBRS raw files are 5-13M rows per year. `tidyverse` is acceptable for plotting and small post-estimation work.
- **Identifiers**: 5-digit county FIPS as character (zero-padded). State FIPS as 2-digit character.
- **Standard errors**: Clustered at county level by default. Two-way (county + year) only when explicitly requested. Type HC1.
- **Significance stars**: `*` p<0.10, `**` p<0.05, `***` p<0.01.

## Directory layout

```
green-borders/
  CLAUDE.md                       # this file
  README.md                       # human-facing overview
  .gitignore                      # data/raw is gitignored
  renv.lock                       # pinned package versions
  R/
    00_master.R                   # sources 00 -> 09 in order
    00_build_metadata.R           # county metadata, treated assignment
    01_nibrs_load.R               # NIBRS Kaplan V9 -> county-year panel
    01c_iowadot_load.R            # Iowa DOT PDFs -> county-year panel
    02_acs.R                      # Census API -> demographics
    03_panel.R                    # merge all sources -> analytical panel
    04_raw_plot.R                 # Figure 2 raw trends
    04_gravity_regression.R       # Chapter 4 PPML gravity (TODO)
    05_twfe.R                     # baseline DiD
    06_eventstudy.R               # event study by outcome
    07_robust.R                   # robustness specs
    08_iowadot_did.R              # DiD on Iowa DOT outcomes
    09_sumstats.R                 # summary stats, balance table
    utils/
      twfe_no_fixest.R            # lm + sandwich fallback
      spatial_aggregation.R       # radius population helpers (TODO)
      pdf_parser.R                # Iowa DOT PDF -> tidy long
  data/
    raw/                          # GITIGNORED; original downloads
      nibrs/                      # Kaplan V9 per-year .rds
      leaic/                      # ORI -> county FIPS crosswalk .rda
      iowadot/                    # Iowa DOT PDFs
      shapefiles/                 # county polygons
    interim/                      # filtered, normalized intermediates
    processed/                    # final analytical panel + metadata
  output/
    tables/                       # .tex, .csv
    figures/                      # .pdf, .png
  docs/
    research_design.md            # specification of clusters 1-11
    feedback_lee.md               # Tom Lee defense feedback (6 items)
    feedback_lychagin.md          # Sergey meeting reconstruction
  tests/                          # unit tests for utils/
  paper/                          # thesis manuscript (LaTeX)
```

## Research design

The empirical specifications are organized into numbered "clusters" in `docs/research_design.md`. When implementing a feature, **reference cluster numbers** so the lineage is explicit. The clusters are:

- **Cluster 1**: Treatment definition. IL-border counties in IA, IN, WI as treated; same-state interior as control; MO as placebo (Feb 2023 own legalization).
- **Cluster 2**: Sample window 2015-2022 (NIBRS Kaplan V9 limit). Phase 2 will extend to 2023-2024 once V11 is available.
- **Cluster 4**: Outcomes. Drug arrests (headline), OWI/DUI (NIBRS Group B), property crime, violent crime. Iowa DOT OWI revocations + speeding convictions as registry-based corroboration.
- **Cluster 5**: Baseline TWFE specification with county FE + state-by-year FE.
- **Cluster 6**: Gravity entry regression for Chapter 4. PPML with three radii (15, 50, 100 miles).
- **Cluster 7**: Event study with k = -5..2 around 2020, k = -1 as reference.
- **Cluster 8**: Robustness battery (drop 2020 COVID, drop MO, drop largest treated counties, population-weighted, per-state).
- **Cluster 10**: Missouri 2023 placebo as secondary identification check.
- **Cluster 11**: Cross-border population interaction terms for the gravity model.

## Tasks Claude Code may be given

Tasks should always reference (a) a cluster number, (b) the script file(s) to read or modify, and (c) explicit input/output expectations. Examples are in `docs/claude_code_task_examples.md`.

## What Claude Code should not do

- Do not commit anything in `data/raw/`. Verify before any `git add` that nothing under that directory is staged.
- Do not change the cluster definitions in `docs/research_design.md` without an explicit instruction from the user.
- Do not auto-bump R package versions in `renv.lock`. The environment is fragile (rlang 1.1.6 vs 1.1.7 incompatibility with arrow/fixest/tidycensus); package changes go through a separate workflow.
- Do not use `arrow::write_parquet` or `arrow::read_parquet`. Use `saveRDS` / `readRDS` with `.rds` extension. The arrow package has the same rlang incompatibility.
- Do not bring in the Census API key as a literal. Use `Sys.getenv("CENSUS_API_KEY")` and assume the user has set it in `.Renviron`.

## Style notes for code generation

- Use `here::here()` for all file paths.
- Function names: snake_case verb-first (`compute_radius_population`, `parse_iadot_pdf`).
- Defensive about column existence; use the `pick_col(dt, candidates)` pattern from `R/01_nibrs_load.R` when a column might appear under several names.
- Include `cat()` progress prints in long-running scripts (NIBRS load, gravity).
- Every analysis script (`04_*` through `09_*`) writes both a `.csv` and a `.tex` for every table.

## Known data quirks

- **NIBRS Kaplan V9** uses lowercase state values (`"iowa"`, not `"IOWA"`), English descriptions for `ucr_offense_code` (`"drug/narcotic offenses - drug/narcotic violations"`, not the NIBRS code `"35A"`), and stores OWI under `ucr_arrest_offense_code` in the Group B segment, not in the offense segment.
- **2018 NIBRS coverage gap**: Iowa-only event-study showed large negative outliers at k=-2 (2018) across all four outcomes. This is a reporting agency coverage rotation artifact, not a behavioral pre-trend. Note it in any robustness commentary.
- **Iowa DOT PDFs**: county names are in uppercase; "O'Brien" appears as "OBRIEN" (no apostrophe).
