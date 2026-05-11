# Advisor Feedback: Sergey Lychagin (Meeting Reconstruction)

Implementation-focused summary of advisor meeting. Discussion of advising relationship and synthesis are not reproduced here; this document is what Claude Code reads.

## Repositioning

- Crime is the **primary outcome** (job-market paper). Dispensary entry is supporting / Chapter 4.
- Narrative arc: entry (lead-up) → crime (main) → federalism externalities (bigger story).

## DiD setup

- Treated: counties bordering Illinois in prohibition states (IA, IN, KY, WI). KY excluded in practice due to data sparsity (see `docs/research_design.md` cluster 1).
- Control: same-state interior counties.
- Three candidate control groups to report separately:
  1. Interior counties in the same prohibition state.
  2. IL counties bordering legal states (WI, MI).
  3. Counties in non-IL prohibition states bordering other prohibition states (placebo).
- "Separate regressions = every regression": each control group gets its own column.

## Gravity entry regression

- PPML, three radii: 15-mile (interaction term), 50-mile (sensitivity), 100-mile (gravity baseline).
- Specification:
  `Entry = exp(α + β·Border_100 + log(Pop_IL) + log(Pop_PROH) + ε)`
- Pop_IL and Pop_PROH must be measured **separately**; do not collapse into a single "neighbor pop" variable.
- Missouri 2023 held out from primary identification.

## Crime DiD

- Outcomes: DUI / OWI, drug possession, vehicular incidents.
- Specification: county FE + state-by-year FE (not just year FE) to absorb each prohibition state's idiosyncratic trajectory.
- Cluster SE at county level.
- Report event study with pre-trends explicitly.
- Iowa border counties drive identification; expansion to IN + WI in Phase 2 boosts power.

## Phase 2 (post-defense)

- NIBRS V11 covers 2023-2024. Re-estimate with Missouri 2023 as a second event.
- FOIA Iowa State Patrol for cannabis-specific interdiction statistics by district-year.
- Consider extending to Michigan 2018 as an earlier shock.

## Top quotable

> "Fiscal forecasts must net out borrowed rents."

(Candidate one-liner for thesis abstract and policy slide.)
