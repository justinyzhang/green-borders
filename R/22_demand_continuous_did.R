# =============================================================================
# 22_demand_continuous_did.R
# =============================================================================
# Chapter 5 demand-side analysis. Replaces the v15 binary treatment with a
# CONTINUOUS treatment intensity equal to Illinois-side population within
# radius r of each Iowa county centroid. Three specifications:
#
#   (A) Continuous-dose DiD a la Callaway-Goodman-Bacon-Sant'Anna (NBER WP
#       32117). Treatment = log(PopIL_{c,r} + 1) varies across all 99 IA
#       counties; identification leverages within-county variation in dose
#       interacted with the post-2020 indicator.
#
#   (B) Holmes-style banded DiD. Partition IA counties into bands by their
#       Illinois-side catchment population at three radii (15, 50, 100 mi),
#       interact each band with Post. Prediction: monotone decrease in
#       coefficient as band moves outward, with the bulk of the effect in
#       the innermost band.
#
#   (C) Continuous-dose event study. Treatment x year interactions at each
#       k = -5..+2 relative to 2020.
#
# Comparison: the v15 binary headline yields beta = 0.422 (SE 0.198,
# p = 0.033) for log(drug arrests + 1). When dose is evaluated at the median
# treated bridge-county dose, the continuous spec should recover something
# close to this; the continuous spec provides much tighter identification
# because it uses 99 county clusters of variation rather than 8 vs. 91.
#
# Inputs:
#   - data/iowa_panel.csv (existing: 99 IA counties x 8 years, NIBRS arrests
#     by category + controls). Columns required:
#       GEOID, county, year, drug_arrests, property_arrests, violent_arrests,
#       owi_arrests, population, income
#   - out/catchment_pop_ia_focal.csv (from 20_catchment_populations.R)
#
# Outputs:
#   - out/demand_continuous_results.csv
#   - out/demand_banded_results.csv
#   - out/demand_eventstudy_results.csv
#   - out/demand_table.tex (booktabs LaTeX for the paper)
#   - out/fig4_continuous_eventstudy.pdf (AER-style event study plot)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
})

OUT_DIR <- "out"
dir.create(OUT_DIR, showWarnings = FALSE)

# AER plot palette (matches fig3_eventstudy_iowa.pdf in v15 appendix)
AER_RED  <- "#9E2A2B"
AER_GREY <- "#666666"

# ---- 1. Load Iowa panel and catchments -------------------------------------

ia_panel <- read_csv("data/iowa_panel.csv", show_col_types = FALSE)
stopifnot(nrow(ia_panel) == 99 * 8, all(c("GEOID","year","drug_arrests") %in% names(ia_panel)))

cat_ia <- read_csv("out/catchment_pop_ia_focal.csv", show_col_types = FALSE)

# Reshape catchments to wide: one row per county with pop_IL at each radius
cat_ia_wide <- cat_ia %>%
  pivot_wider(
    id_cols = c(GEOID, COUNAME),
    names_from = radius_mi,
    values_from = c(pop_IL, pop_IA),
    names_glue = "{.value}_{radius_mi}mi"
  )

# Merge into panel
ia <- ia_panel %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
  left_join(cat_ia_wide, by = "GEOID") %>%
  mutate(
    post = as.integer(year >= 2020),
    log_drug    = log(drug_arrests    + 1),
    log_prop    = log(property_arrests + 1),
    log_viol    = log(violent_arrests  + 1),
    log_owi     = log(owi_arrests      + 1),
    log_pop     = log(population),
    log_inc     = log(income),
    # Continuous treatment intensity at each radius (log Illinois-side pop)
    dose_15  = log(pop_IL_15mi  + 1),
    dose_50  = log(pop_IL_50mi  + 1),
    dose_100 = log(pop_IL_100mi + 1),
    # Center doses at the within-Iowa mean for interpretation. The level
    # coefficient on dose:post is then the differential effect of a one-unit
    # increase in log-IL-pop relative to the average Iowa county.
    dose_15_c  = dose_15  - mean(dose_15,  na.rm = TRUE),
    dose_50_c  = dose_50  - mean(dose_50,  na.rm = TRUE),
    dose_100_c = dose_100 - mean(dose_100, na.rm = TRUE)
  )

# Sanity: confirm bridge counties have highest doses
bridge_geoids <- c("19045","19057","19061","19097","19111","19139","19163","19115")
cat("\nDose values for v15 treated (bridge) counties at r = 50 mi:\n")
print(ia %>%
        filter(GEOID %in% bridge_geoids, year == 2019) %>%
        select(GEOID, county, pop_IL_50mi, dose_50) %>%
        arrange(desc(pop_IL_50mi)))

cat("\nDose summary across all 99 IA counties at r = 50 mi:\n")
print(summary(ia$dose_50[ia$year == 2019]))

# ---- 2. Specification A: Continuous-dose DiD -------------------------------

# For each radius, run: outcome = alpha_c + lambda_t + delta * dose_c * post_t + X_ct b
# Cluster SE at county level (99 clusters; CGS asymptotics now well-justified).

cat("\n=== Specification A: Continuous-dose DiD (drug arrests, 4 outcomes) ===\n")

run_cont_did <- function(outcome, dose_var, data) {
  fml <- as.formula(sprintf(
    "%s ~ i(post, %s, ref = 0) + log_pop + log_inc | GEOID + year",
    outcome, dose_var
  ))
  feols(fml, data = data, cluster = ~GEOID, lean = TRUE)
}

# Run for the three radii and four outcomes
outcomes <- c("log_drug","log_prop","log_viol","log_owi")
radii    <- c("dose_15","dose_50","dose_100")

cont_models <- list()
for (out in outcomes) {
  for (d in radii) {
    label <- paste0(out, "__", d)
    cont_models[[label]] <- run_cont_did(out, d, ia)
  }
}

# Extract dose:post coefficient from each model
extract_cont <- function(m, label) {
  s <- summary(m)
  ct <- s$coeftable
  # Find the dose:post interaction row (fixest names it like "post::1:dose_50")
  i <- grep("post::1:dose", rownames(ct))
  if (length(i) == 0) return(NULL)
  tibble(
    spec = label,
    beta = ct[i, "Estimate"],
    se   = ct[i, "Std. Error"],
    t    = ct[i, "t value"],
    p    = ct[i, "Pr(>|t|)"],
    N    = s$nobs,
    n_clusters = length(unique(s$fixef_id[[1]]))
  )
}
cont_results <- imap_dfr(cont_models, ~ extract_cont(.x, .y)) %>%
  separate(spec, into = c("outcome","radius"), sep = "__") %>%
  mutate(
    outcome_label = recode(outcome,
                           log_drug = "Drug",   log_prop = "Property",
                           log_viol = "Violent", log_owi  = "OWI"),
    radius_label  = recode(radius,
                           dose_15 = "15 mi", dose_50 = "50 mi", dose_100 = "100 mi")
  )

cat("\nContinuous-dose DiD coefficients (dose:post interaction):\n")
print(cont_results %>%
        select(outcome_label, radius_label, beta, se, p, N) %>%
        mutate(across(c(beta, se), ~ round(., 4)),
               p = round(p, 3)))

write_csv(cont_results, file.path(OUT_DIR, "demand_continuous_results.csv"))

# ---- 3. Sanity check: compare to v15 binary specification ------------------

# At dose = log(median treated bridge county's IL-pop within 50mi),
# the continuous specification should imply approximately the v15 binary 0.422.

cat("\n=== Sanity check: continuous-dose vs. v15 binary ===\n")

# Run v15 binary spec on the same panel
ia <- ia %>%
  mutate(treated_v15 = as.integer(GEOID %in% bridge_geoids))

v15_check <- feols(log_drug ~ i(post, treated_v15, ref = 0) + log_pop + log_inc |
                     GEOID + year,
                   data = ia, cluster = ~GEOID, lean = TRUE)
cat("v15 binary headline replication: beta =",
    round(coef(v15_check)["post::1:treated_v15"], 4),
    " SE =", round(se(v15_check)["post::1:treated_v15"], 4), "\n")
cat("(Compare to v15 reported: beta = 0.422, SE = 0.198)\n")

# Implied effect at median bridge dose using continuous spec at r = 50 mi
bridge_dose_median <- ia %>%
  filter(GEOID %in% bridge_geoids, year == 2019) %>%
  summarise(med = median(dose_50, na.rm = TRUE)) %>%
  pull(med)

interior_dose_median <- ia %>%
  filter(!(GEOID %in% bridge_geoids), year == 2019) %>%
  summarise(med = median(dose_50, na.rm = TRUE)) %>%
  pull(med)

dose_50_delta <- bridge_dose_median - interior_dose_median
cont_drug_50_beta <- cont_results %>%
  filter(outcome == "log_drug", radius == "dose_50") %>%
  pull(beta)
implied_effect <- cont_drug_50_beta * dose_50_delta
cat("Continuous-dose implied differential (bridge median - interior median) at r=50:\n")
cat("  dose differential = log-IL-pop, treated median - control median =",
    round(dose_50_delta, 3), "\n")
cat("  beta_continuous * delta = ", round(cont_drug_50_beta, 4), " *",
    round(dose_50_delta, 3), " =",
    round(implied_effect, 4), "\n")
cat("(Should be in the same neighborhood as v15 binary 0.422.)\n")

# ---- 4. Specification B: Holmes-style banded DiD ---------------------------

# Partition IA counties into bands by their IL-side catchment population at r=50.
# Band 1 (highest): top quintile of pop_IL_50mi (innermost 20%)
# Band 2: 2nd quintile
# Band 3: 3rd quintile
# Band 4: 4th quintile
# Band 5 (reference): bottom quintile (essentially zero IL exposure)

cat("\n=== Specification B: Holmes-style banded DiD ===\n")

ia <- ia %>%
  group_by(year) %>%
  mutate(band_50 = ntile(pop_IL_50mi, 5)) %>%
  ungroup() %>%
  mutate(band_50 = factor(band_50, levels = 1:5,
                          labels = c("Q1 (lowest)","Q2","Q3","Q4","Q5 (highest)")))

# Verify: bridge counties should be in the top band
cat("\nBand distribution for v15 bridge counties:\n")
print(table(ia$band_50[ia$year == 2019 & ia$GEOID %in% bridge_geoids]))

# Run banded DiD with band 1 (Q1 lowest) as reference
banded_drug <- feols(
  log_drug ~ i(band_50, post, ref = "Q1 (lowest)") + log_pop + log_inc |
    GEOID + year,
  data = ia, cluster = ~GEOID, lean = TRUE
)
summary(banded_drug)

# Extract banded coefficients
banded_results <- broom::tidy(banded_drug, conf.int = TRUE) %>%
  filter(grepl("band_50.*post", term)) %>%
  mutate(band = gsub(".*::(Q[1-5][^:]*).*", "\\1", term))

cat("\nBanded DiD coefficients (drug arrests):\n")
print(banded_results %>% select(band, estimate, std.error, p.value))

write_csv(banded_results, file.path(OUT_DIR, "demand_banded_results.csv"))

# ---- 5. Specification C: Continuous-dose event study -----------------------

cat("\n=== Specification C: Continuous-dose event study ===\n")

# Event time relative to 2020. Reference: k = -1 (2019).
ia <- ia %>%
  mutate(event_time = year - 2020)

run_event_study <- function(outcome, dose_var, data, label) {
  fml <- as.formula(sprintf(
    "%s ~ i(event_time, %s, ref = -1) + log_pop + log_inc | GEOID + year",
    outcome, dose_var
  ))
  m <- feols(fml, data = data, cluster = ~GEOID, lean = TRUE)
  # Tidy: extract event-time coefficients
  tidy_m <- broom::tidy(m, conf.int = TRUE, conf.level = 0.95) %>%
    filter(grepl("event_time::", term)) %>%
    mutate(
      k = as.integer(gsub("event_time::(-?[0-9]+):.*", "\\1", term)),
      outcome = label
    ) %>%
    bind_rows(tibble(k = -1, estimate = 0, std.error = 0,
                     conf.low = 0, conf.high = 0, outcome = label)) %>%
    arrange(k)
  list(model = m, tidy = tidy_m)
}

es_drug_50 <- run_event_study("log_drug", "dose_50", ia, "Drug arrests")
es_drug_15 <- run_event_study("log_drug", "dose_15", ia, "Drug arrests (r=15)")
es_drug_100 <- run_event_study("log_drug", "dose_100", ia, "Drug arrests (r=100)")

cat("\nDrug arrests event study (r = 50 mi, dose-scaled coefficients):\n")
print(es_drug_50$tidy %>% select(k, estimate, std.error, conf.low, conf.high) %>%
        mutate(across(estimate:conf.high, ~ round(., 4))))

# Joint Wald test of pre-period coefficients
pre_coefs_idx <- grep("event_time::-[2-5]:", names(coef(es_drug_50$model)))
if (length(pre_coefs_idx) > 0) {
  pre_coef_names <- names(coef(es_drug_50$model))[pre_coefs_idx]
  wald_test <- wald(es_drug_50$model, keep = pre_coef_names, print = FALSE)
  cat("\nJoint Wald test of pre-period coefficients (k = -5, -4, -3, -2):\n")
  cat("  F =", round(wald_test$stat, 3),
      "  p_joint =", round(wald_test$p, 4), "\n")
}

# Combine all event study results
es_all <- bind_rows(
  es_drug_15$tidy %>% mutate(radius = "15 mi"),
  es_drug_50$tidy %>% mutate(radius = "50 mi"),
  es_drug_100$tidy %>% mutate(radius = "100 mi")
)
write_csv(es_all, file.path(OUT_DIR, "demand_eventstudy_results.csv"))

# ---- 6. AER-style event study plot -----------------------------------------

# Use the same conventions as fig3_eventstudy_iowa.pdf:
# - red points & connecting line
# - grey 95% CI error bars
# - solid y=0 reference
# - dotted x=0 reference
# - light grey gridlines
# - in-plot panel label

es_plot_data <- es_drug_50$tidy

p_es <- ggplot(es_plot_data, aes(x = k, y = estimate)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted", linewidth = 0.4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.15, color = AER_GREY, linewidth = 0.5) +
  geom_line(color = AER_RED, linewidth = 0.6) +
  geom_point(color = AER_RED, size = 2.4) +
  scale_x_continuous(breaks = -5:2) +
  labs(
    x = "Years from January 2020 IL legalization",
    y = "Dose-scaled coefficient on log(IL pop within 50 mi)",
    title = NULL
  ) +
  annotate("text", x = -5, y = max(es_plot_data$conf.high) * 0.95,
           label = "Panel A. Continuous-dose drug arrests (r = 50 mi)",
           hjust = 0, size = 4) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

ggsave(file.path(OUT_DIR, "fig4_continuous_eventstudy.pdf"),
       p_es, width = 7, height = 4.2)
cat("\nWrote", file.path(OUT_DIR, "fig4_continuous_eventstudy.pdf"), "\n")

# ---- 7. LaTeX table assembly -----------------------------------------------

# Main table: 4 outcomes x 3 radii continuous-dose results, headline column
# being drug arrests at r = 50 mi.

main_models <- list(
  "Drug, r=15"   = cont_models[["log_drug__dose_15"]],
  "Drug, r=50"   = cont_models[["log_drug__dose_50"]],
  "Drug, r=100"  = cont_models[["log_drug__dose_100"]],
  "Property, r=50"  = cont_models[["log_prop__dose_50"]],
  "Violent, r=50"   = cont_models[["log_viol__dose_50"]],
  "OWI, r=50"       = cont_models[["log_owi__dose_50"]]
)

ms_demand <- modelsummary(
  main_models,
  output = "latex",
  stars = c("*" = 0.10, "**" = 0.05, "***" = 0.01),
  coef_omit = "log_pop|log_inc",
  coef_rename = c("post::1:dose_15"  = "Dose (15 mi) x Post",
                  "post::1:dose_50"  = "Dose (50 mi) x Post",
                  "post::1:dose_100" = "Dose (100 mi) x Post"),
  gof_omit = "IC|Log|RMSE|Pseudo|Adj.|R2 Within|FE",
  fmt = 3,
  add_rows = data.frame(
    a = c("County FE", "Year FE", "Time-varying controls"),
    b = c("Yes","Yes","Yes"), c = c("Yes","Yes","Yes"),
    d = c("Yes","Yes","Yes"), e = c("Yes","Yes","Yes"),
    f = c("Yes","Yes","Yes"), g = c("Yes","Yes","Yes")
  ),
  notes = paste(
    "Notes: Continuous-dose DiD coefficients on the dose-by-post interaction.",
    "Dose = log(Illinois-side population within r miles of county centroid + 1).",
    "Iowa panel 2015-2022, 99 counties x 8 years = 792 observations.",
    "County and year fixed effects, log population and log income controls.",
    "Standard errors clustered at county level (99 clusters).",
    "Stars: ***p<0.01, **p<0.05, *p<0.10."
  )
)
writeLines(ms_demand, file.path(OUT_DIR, "demand_table.tex"))

cat("\n[DONE]\n")
cat("Outputs in", OUT_DIR, ":\n")
cat("  demand_continuous_results.csv\n")
cat("  demand_banded_results.csv\n")
cat("  demand_eventstudy_results.csv\n")
cat("  fig4_continuous_eventstudy.pdf\n")
cat("  demand_table.tex\n")
