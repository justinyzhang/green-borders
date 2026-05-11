# Research Design

This document specifies the empirical implementation of the project. Specifications are grouped into numbered clusters. Claude Code tasks reference these cluster numbers.

## Cluster 1 — Treatment definition

**Treated units**: Counties in prohibition states that share any land or river border with the state of Illinois.

**Sample states (3)**: Iowa (IA), Indiana (IN), Wisconsin (WI).

**Excluded states**:
- **Missouri (MO)**: legalized recreational cannabis Feb 2023; held out from main panel and used as a secondary placebo (cluster 10).
- **Kentucky (KY)**: only Hancock County borders Illinois (across the Ohio River), and KY NIBRS coverage pre-2020 is sparse. Excluded.
- **Michigan (MI)**: legalized 2018; treated as a separate legalization shock, not part of this design.

**Treated county FIPS** (verified from US Census Tiger):

| State | FIPS | Counties |
|-------|------|----------|
| IA (10) | 19005, 19043, 19045, 19057, 19061, 19097, 19111, 19115, 19139, 19163 | Allamakee, Clayton, Des Moines, Dubuque, Henry, Jackson, Lee, Louisa, Muscatine, Scott |
| WI (5) | 55045, 55059, 55101, 55105, 55127 | Green, Kenosha, Racine, Rock, Walworth |
| IN (12) | 18007, 18011, 18051, 18083, 18091, 18127, 18129, 18167, 18171, 18175, 18179, 18181 | Benton, Boone-area, Gibson, Knox, Lake, Newton, Porter, Vermillion, Vigo, Warren, Warrick, White |

Total treated: **27 counties**.

**Control units**: Same-state non-IL-border counties in IA + IN + WI. Total: **236 counties** (89 IA + 67 WI + 80 IN). Note: WI has 72 counties total, IN has 92, IA has 99.

**Total main panel**: 263 counties × 8 years = 2,104 obs.

**Missouri placebo set**: 10 IL-bordering counties (29007, 29019, 29045, 29113, 29127, 29163, 29173, 29183, 29186, 29199, 29219, 29510) + 105 MO interior. Used only in cluster 10.

## Cluster 2 — Time window

**Main sample**: 2015-2022 (8 years). Bounded above by Kaplan NIBRS V9 release (ends 2022).

**Treatment timing**: Post = 1 if `year >= 2020` (Illinois sales began Jan 1, 2020).

**Pre-period**: 2015-2019 (5 years).
**Post-period**: 2020-2022 (3 years).

**Phase 2 extension** (not in defense draft): When NIBRS V11 becomes available (covers 2023-2024), extend the panel and add the Missouri 2023 shock as a staggered second event.

## Cluster 4 — Outcomes

Four NIBRS-based outcomes, all aggregated to county-year level.

| Outcome | Source segment | Column in panel |
|---------|----------------|-----------------|
| Drug arrests (headline) | offense_segment | `drug` |
| OWI / DUI arrests | group_b_arrest_report_segment | `owi` |
| Property crime | offense_segment | `property` |
| Violent crime | offense_segment | `violent` |

Each outcome enters the regression as `ln(Y + 1)`. Rates are also computed as `Y / pop_total * 1e5` and stored in `*_rate` columns for plotting.

**Drug** = sum of `"drug/narcotic offenses - drug/narcotic violations"` and `"drug/narcotic offenses - drug equipment violations"`.

**OWI** = sum across `"driving under the influence"` and synonyms in the Group B `ucr_arrest_offense_code` field. See `R/01_nibrs_load.R` for the full synonym list.

**Property** = burglary, larceny (all subtypes), motor vehicle theft, destruction/damage of property.

**Violent** = aggravated/simple/intimidation assault, robbery, murder/manslaughter, rape and other sex offenses, kidnapping.

### Iowa DOT corroboration outcomes (cluster 4b)

When Iowa DOT registry data is available, additional outcomes:

| Outcome | Source | File |
|---------|--------|------|
| OWI revocations | Iowa DOT MVD (PDF, 2005-2024) | `data/raw/iowadot/owi_revocations.pdf` |
| Speeding convictions | Iowa DOT MVD (PDF, 2013-2022) | `data/raw/iowadot/speeding_convictions.pdf` |
| Traffic fatalities | Iowa DOT (PDF, 2015-2024) | `data/raw/iowadot/yearly_fatalities.pdf` |

Iowa DOT data is Iowa-only and serves as a registry-based cross-check on the NIBRS results. Useful because it bypasses NIBRS reporting agency coverage rotation.

## Cluster 5 — Baseline specification

```
ln(Y_cst + 1) = alpha_c + lambda_st + tau * (Border_c * Post_t)
              + gamma_1 * ln_pop_cst + gamma_2 * ln_inc_cst
              + u_cst
```

where `c` indexes county, `s` state, `t` year.

- `alpha_c`: county fixed effect.
- `lambda_st`: **state-by-year fixed effect** (not just year FE; absorbs each state's idiosyncratic time path).
- `Border_c`: 1 if county c is IL-border (cluster 1 list).
- `Post_t`: 1 if `year >= 2020`.
- `ln_pop`, `ln_inc`: log ACS population, log median household income.
- Standard errors clustered at county level. Type HC1.

Also report a **Spec B** with `lambda_t` (year FE only) for comparison.

The treatment indicator `treat = Border * Post` is precomputed in `03_panel.R`.

## Cluster 6 — Gravity entry regression (Chapter 4)

PPML specification for IL dispensary entry, replacing the OLS specification in defense v4 slides.

```
Entry_i = exp(beta_0 + beta_1 * Border_15_i + beta_2 * log(Pop_IL_i)
               + beta_3 * log(Pop_PROH_i) + epsilon_i)
```

- Unit: 102 Illinois counties.
- `Entry_i`: count of operational dispensaries by some date (e.g., June 2024).
- `Border_15_i`: 1 if county i has any boundary point within 15 miles of a non-legal state.
- `Pop_IL_i`: population of IL within X miles of i (Census 2020).
- `Pop_PROH_i`: population of prohibition states within X miles of i.
- Estimated three times at X = 15, 50, 100 miles (sensitivity).
- Clustered SE at county.
- Estimator: `fixest::fepois`, or `glm(family = quasipoisson)` if fixest unavailable.

Report three columns in the main table:
1. Binary border dummy only (baseline).
2. Border dummy × 15-mile prohibition population (interaction).
3. 15-mile prohibition population in levels.

## Cluster 7 — Event study

```
ln(Y_cst + 1) = sum_{k != -1} beta_k * 1{year_rel == k} * Border_c
              + alpha_c + lambda_st + gamma * ln_pop + u_cst
```

- `year_rel = year - 2020`.
- Event window: `k = -5, -4, -3, -2, -1, 0, 1, 2`.
- Reference period: `k = -1` (2019).
- One event study per outcome (drug, OWI, property, violent).
- 95% CIs from cluster-robust SE.
- Plot using `ggplot2` with `geom_errorbar` + `geom_point`.

## Cluster 8 — Robustness

Apply to the **drug arrests** headline regression. Each spec is a row in `output/tables/tab2_robust.tex`.

1. **Baseline**: cluster 5 reproduction.
2. **Drop 2020**: remove year 2020 entirely (COVID donut).
3. **Drop Lake County, IN (18089)**: Chicago suburb confound.
4. **Drop Scott County, IA (19163)**: 35% of IA treated population.
5. **Drop Missouri-border-via-Mississippi**: keep only WI + IN land-border treated counties.
6. **Population-weighted**: `weights = pop_total`.
7. **State-by-state**: separate IA-only, IN-only, WI-only regressions.

For each spec, report β on `treat`, cluster-robust SE, p-value, N obs, N counties.

## Cluster 10 — Missouri 2023 placebo

Secondary identification check. Available only when NIBRS V11 is in hand (Phase 2). Implementation stub:

```
ln(Y_cst + 1) = alpha_c + lambda_t
              + tau_MO * (MO_Border_c * Post_2023_t)
              + u_cst
```

- Sample: MO counties (115 total). 10 IL-border + 105 interior.
- `Post_2023_t = 1` if `year >= 2023`.
- If `tau_MO > 0` significantly, this strengthens the cluster-5 result by showing the same shock (now MO->KS/AR/OK/TN) produces a similar pattern.

## Cluster 11 — Cross-border population interaction

Extension to cluster 5. Replaces binary `Border_c * Post_t` with continuous interaction:

```
ln(Y_cst + 1) = alpha_c + lambda_st + tau * (PopIL_15mi_c * Post_t)
              + gamma_1 * ln_pop + gamma_2 * ln_inc + u_cst
```

- `PopIL_15mi_c`: log IL population within 15 miles of county c centroid.
- This sharpens the dose-response interpretation: counties closer to denser IL settlements should show larger treatment effects if the mechanism is cross-border travel.

Helper function in `R/utils/spatial_aggregation.R`:

```r
compute_radius_population(county_fips, radius_miles, state_subset)
```

- `county_fips`: character, 5-digit FIPS.
- `radius_miles`: numeric.
- `state_subset`: character vector of state FIPS codes to sum across (e.g., `c("17")` for IL only).
- Returns: numeric, total population from counties whose centroids fall within `radius_miles` of the given county's centroid.
