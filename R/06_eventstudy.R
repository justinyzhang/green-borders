# =============================================================================
# 06_eventstudy.R   (4-state, state-by-year FE, all 4 outcomes)
#
# Specification:
#   ln(Y_cst + 1) = sum_{k != -1} beta_k * 1{year_rel == k} * Border_c
#                + alpha_c + lambda_st + gamma * ln_pop_cst + u_cst
#
# Sample: 263 counties in IA, IN, WI (Missouri excluded).
#
# Outputs:
#   out/fig3_eventstudy.pdf            (drug, headline)
#   out/figA2_eventstudy_owi.pdf
#   out/figA3_eventstudy_property.pdf
#   out/figA4_eventstudy_violent.pdf
#   out/fig3_eventstudy_coefs.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(here)
})

panel <- readRDS(here::here("data", "processed", "panel_county_year.rds"))
setDT(panel)
panel_main <- panel[treatment %in% c("IL_border", "interior")]
panel_main[, year_rel := year - 2020]
event_years <- sort(unique(panel_main$year_rel))
ref_k <- -1

fit_event <- function(y_col) {
  dt <- copy(panel_main)
  for (k in event_years) {
    if (k == ref_k) next
    col <- sprintf("k_%s", ifelse(k < 0, paste0("m", abs(k)), as.character(k)))
    dt[, (col) := as.integer(year_rel == k) * is_il_border]
  }
  k_cols <- grep("^k_", names(dt), value = TRUE)
  f <- as.formula(sprintf(
    "%s ~ %s + ln_pop + factor(county_fips) + factor(state_year)",
    y_col, paste(k_cols, collapse = " + ")
  ))
  m <- lm(f, data = dt)
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  rows <- data.table(
    k = sapply(k_cols, function(s) {
      x <- sub("^k_", "", s)
      if (substr(x, 1, 1) == "m") -as.integer(sub("^m", "", x)) else as.integer(x)
    }),
    est   = ct[k_cols, "Estimate"],
    se    = ct[k_cols, "Std. Error"],
    pval  = ct[k_cols, "Pr(>|t|)"],
    ci_lo = ci[k_cols, 1],
    ci_hi = ci[k_cols, 2]
  )
  rows <- rbindlist(list(rows, data.table(k = ref_k, est = 0, se = 0,
                                          pval = NA_real_, ci_lo = 0, ci_hi = 0)))
  setorder(rows, k)
  rows
}

plot_es <- function(es, title, ylab, color = "#1f78b4") {
  ggplot(es, aes(k, est)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "gray60") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15,
                  color = color, linewidth = 0.6) +
    geom_point(size = 2.5, color = color) +
    scale_x_continuous(breaks = event_years) +
    labs(x = "Years relative to January 2020 IL legalization",
         y = ylab, title = title,
         subtitle = "Sample: IA + IN + WI; reference period 2019 (k = -1)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)

# Headline: drug
es_drug <- fit_event("ln_drug")
cat("=== Event study: DRUG (headline) ===\n"); print(es_drug)
ggsave(here::here("out", "fig3_eventstudy.pdf"),
       plot_es(es_drug,
               "Event study: drug arrests, IL-bordering vs interior counties (IA+IN+WI)",
               "Coefficient on log drug arrests (95% CI)", "#1f78b4"),
       width = 7.5, height = 4.5)

# Appendix figures
es_owi      <- fit_event("ln_owi")
es_property <- fit_event("ln_property")
es_violent  <- fit_event("ln_violent")

cat("\n=== Event study: OWI ===\n"); print(es_owi)
cat("\n=== Event study: PROPERTY ===\n"); print(es_property)
cat("\n=== Event study: VIOLENT ===\n"); print(es_violent)

ggsave(here::here("out", "figA2_eventstudy_owi.pdf"),
       plot_es(es_owi, "Event study: OWI arrests",
               "Coefficient on log OWI arrests (95% CI)", "#7a7a7a"),
       width = 7.5, height = 4.5)
ggsave(here::here("out", "figA3_eventstudy_property.pdf"),
       plot_es(es_property, "Event study: property crime",
               "Coefficient on log property crime (95% CI)", "#7a7a7a"),
       width = 7.5, height = 4.5)
ggsave(here::here("out", "figA4_eventstudy_violent.pdf"),
       plot_es(es_violent, "Event study: violent crime",
               "Coefficient on log violent crime (95% CI)", "#7a7a7a"),
       width = 7.5, height = 4.5)

all_es <- rbindlist(list(
  data.table(outcome = "drug",     es_drug),
  data.table(outcome = "owi",      es_owi),
  data.table(outcome = "property", es_property),
  data.table(outcome = "violent",  es_violent)
))
fwrite(all_es, here::here("out", "fig3_eventstudy_coefs.csv"))
cat("\nWrote 4 event-study PDFs + fig3_eventstudy_coefs.csv\n")
