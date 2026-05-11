# green-borders

Cross-border demand and crime spillovers from Illinois cannabis legalization.

This repository contains the empirical pipeline for a job-market paper / MA thesis examining (i) the spatial pattern of Illinois recreational cannabis dispensary entry near prohibition-state borders, and (ii) downstream crime spillovers in IL-bordering counties of Iowa, Indiana, and Wisconsin around the January 2020 IL legalization shock.

## What's here

- **`R/`** — analysis scripts, numbered `00_*` through `09_*`. `R/00_master.R` runs the full pipeline.
- **`docs/research_design.md`** — specification of the design clusters referenced throughout. Read this before modifying any script.
- **`docs/feedback_lee.md`** and **`docs/feedback_lychagin.md`** — committee and advisor feedback driving the current design.
- **`CLAUDE.md`** — instructions for Claude Code agents.
- **`docs/claude_code_task_examples.md`** — examples of well-scoped tasks.

## What's not here

- **Raw NIBRS data** (`data/raw/nibrs/`) — Kaplan V9 per-year `.rds` files, ~600 MB. Download from ICPSR study 38649.
- **Raw LEAIC** (`data/raw/leaic/`) — ORI to county FIPS crosswalk. Download from ICPSR study 35158.
- **Iowa DOT PDFs** (`data/raw/iowadot/`) — download from `https://iowadot.gov/mvd/FactsandStats`.
- **Census API key** — set `CENSUS_API_KEY` in `.Renviron`.

## To reproduce

```r
# Install package set
renv::restore()

# Run pipeline
source("R/00_master.R")
```

## Status

- ✅ Iowa-only NIBRS DiD (Chapter 5 v1)
- 🟡 4-state expansion (IA + IN + WI; MO as placebo)
- 🟡 Iowa DOT registry corroboration (OWI revocations, speeding convictions)
- ❌ Gravity entry regression (Chapter 4 redesign)
- ❌ Missouri 2023 staggered DiD (Phase 2, requires NIBRS V11)

## License

Code: MIT. Data: see source notices.
