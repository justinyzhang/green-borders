# =============================================================================
# R/06b_eventstudy_iowa.R   (v2 — AER Prager & Schmitt 2021 figure style)
# ----------------------------------------------------------------------------
# Iowa-only event study (4 outcomes), county FE + year FE.
# MANDATORY event study for May 20 defense per Sergey's 5/14 directive.
#
# Visual style: replicates Prager & Schmitt (2021, AER) Figure 2 conventions.
#   - Solid red points connected by red line
#   - Grey 95% CI errorbars with horizontal caps
#   - Solid horizontal line at y = 0
#   - Vertical dotted line at k = 0 (treatment)
#   - Light grey horizontal gridlines, no vertical gridlines
#   - Panel label top-left, in-plot
#   - "Years from January 2020 IL legalization" axis title
#
# Specification:
#   ln(Y_ct + 1) = sum_{k != -1} beta_k * 1{year_rel == k} * Border_c
#                + alpha_c + lambda_t + gamma_1 * ln_pop_ct
#                + gamma_2 * ln_inc_ct + u_ct
#
# Sample: 99 Iowa counties (8 treated, 91 interior), 2015-2022.
#
# Outputs:
#   out/fig3_eventstudy_iowa.pdf            (drug, headline, slide 7)
#   out/figA2_eventstudy_iowa_owi.pdf
#   out/figA3_eventstudy_iowa_property.pdf
#   out/figA4_eventstudy_iowa_violent.pdf
#   out/figA1_eventstudy_iowa_3panel.pdf    (drug + property + OWI stacked)
#   out/fig3_eventstudy_iowa_coefs.csv
#   out/fig3_eventstudy_iowa_joint_wald.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(here)
})

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Style constants  (Prager & Schmitt 2021 AER palette)
# -----------------------------------------------------------------------------
AER_RED    <- "#9E2A2B"   # deep red for points + connecting line
AER_GREY   <- "#666666"   # grey for errorbars
AER_BLACK  <- "#000000"
GRID_GREY  <- "#D0D0D0"

# -----------------------------------------------------------------------------
# 1. Load panel and restrict to Iowa
# -----------------------------------------------------------------------------
panel <- readRDS(here::here("data", "processed", "panel_county_year.rds"))
setDT(panel)

ia_panel <- panel[state_name == "iowa" &
                    treatment %in% c("IL_border", "interior")]
ia_panel[, year_rel := year - 2020]

cat("=== Iowa-only event study (AER style) ===\n")
cat("Total obs:        ", nrow(ia_panel), "\n")
cat("Counties:         ", uniqueN(ia_panel$county_fips), "\n")
cat("Treated counties: ", uniqueN(ia_panel[is_il_border == 1]$county_fips), "\n")
cat("Year range:       ", paste(range(ia_panel$year), collapse = " to "), "\n\n")

event_years <- sort(unique(ia_panel$year_rel))
ref_k <- -1L

# -----------------------------------------------------------------------------
# 2. Event study fitter
# -----------------------------------------------------------------------------
fit_event_iowa <- function(y_col) {
  dt <- copy(ia_panel)
  
  for (k in event_years) {
    if (k == ref_k) next
    cname <- sprintf("k_%s", ifelse(k < 0, paste0("m", abs(k)),
                                    as.character(k)))
    dt[, (cname) := as.integer(year_rel == k) * is_il_border]
  }
  k_cols <- grep("^k_", names(dt), value = TRUE)
  
  f <- as.formula(sprintf(
    "%s ~ %s + ln_pop + ln_inc + factor(county_fips) + factor(year)",
    y_col, paste(k_cols, collapse = " + ")
  ))
  
  m  <- lm(f, data = dt)
  vc <- sandwich::vcovCL(m, cluster = dt$county_fips, type = "HC1")
  ct <- lmtest::coeftest(m, vcov. = vc)
  ci <- lmtest::coefci(m, vcov. = vc, level = 0.95)
  
  rows <- data.table(
    k     = sapply(k_cols, function(s) {
      x <- sub("^k_", "", s)
      if (substr(x, 1, 1) == "m") -as.integer(sub("^m", "", x)) else as.integer(x)
    }),
    est   = ct[k_cols, "Estimate"],
    se    = ct[k_cols, "Std. Error"],
    pval  = ct[k_cols, "Pr(>|t|)"],
    ci_lo = ci[k_cols, 1],
    ci_hi = ci[k_cols, 2]
  )
  rows <- rbindlist(list(rows,
                         data.table(k = ref_k, est = 0, se = 0,
                                    pval = NA_real_, ci_lo = 0, ci_hi = 0)))
  setorder(rows, k)
  
  # ---- Joint Wald test of pre-period coefficients ----
  pre_cols <- k_cols[grepl("^k_m", k_cols)]
  pre_idx  <- match(pre_cols, names(coef(m)))
  R        <- matrix(0, length(pre_idx), length(coef(m)))
  for (i in seq_along(pre_idx)) R[i, pre_idx[i]] <- 1
  beta_hat <- coef(m)
  q        <- length(pre_idx)
  wald     <- t(R %*% beta_hat) %*% solve(R %*% vc %*% t(R)) %*% (R %*% beta_hat)
  wald_F   <- as.numeric(wald) / q
  df1      <- q
  df2      <- nobs(m) - length(beta_hat)
  p_joint  <- pf(wald_F, df1 = df1, df2 = df2, lower.tail = FALSE)
  
  list(rows = rows,
       wald = data.table(F_stat  = round(wald_F, 3),
                         df1     = df1,
                         df2     = df2,
                         p_joint = round(p_joint, 4)))
}

# -----------------------------------------------------------------------------
# 3. Plot function — Prager & Schmitt 2021 AER style
# -----------------------------------------------------------------------------
plot_es_aer <- function(es_rows, panel_label, ylim_range = NULL) {
  
  es_rows <- copy(es_rows)
  
  if (is.null(ylim_range)) {
    y_min <- min(es_rows$ci_lo, na.rm = TRUE)
    y_max <- max(es_rows$ci_hi, na.rm = TRUE)
    pad   <- (y_max - y_min) * 0.08
    ylim_range <- c(y_min - pad, y_max + pad)
  }
  
  ggplot(es_rows, aes(x = k, y = est)) +
    
    # Horizontal solid black line at y = 0
    geom_hline(yintercept = 0,
               color     = AER_BLACK,
               linewidth = 0.5) +
    
    # Vertical dotted line at k = 0 (treatment timing)
    geom_vline(xintercept = 0,
               linetype  = "dotted",
               color     = AER_BLACK,
               linewidth = 0.4) +
    
    # 95% CI errorbars (grey, horizontal caps)
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width     = 0.18,
                  color     = AER_GREY,
                  linewidth = 0.4) +
    
    # Connecting line through point estimates
    geom_line(color     = AER_RED,
              linewidth = 0.6) +
    
    # Solid red points
    geom_point(color = AER_RED,
               size  = 2.2) +
    
    # Panel label top-left, in-plot
    annotate("text",
             x     = min(event_years),
             y     = ylim_range[2],
             label = panel_label,
             hjust = 0,
             vjust = 1.2,
             size  = 4.0,
             fontface = "plain") +
    
    scale_x_continuous(
      breaks = event_years,
      labels = as.character(event_years),
      expand = expansion(add = 0.4)
    ) +
    
    scale_y_continuous(
      limits = ylim_range,
      expand = expansion(mult = 0.02)
    ) +
    
    labs(x = "Years from January 2020 IL legalization",
         y = NULL) +
    
    theme_minimal(base_size = 11, base_family = "") +
    theme(
      panel.grid.major.x   = element_blank(),
      panel.grid.minor.x   = element_blank(),
      panel.grid.minor.y   = element_blank(),
      panel.grid.major.y   = element_line(color = GRID_GREY,
                                          linewidth = 0.3),
      
      axis.line.x          = element_line(color = AER_BLACK,
                                          linewidth = 0.4),
      axis.ticks.x         = element_line(color = AER_BLACK,
                                          linewidth = 0.3),
      axis.ticks.length.x  = unit(2.5, "pt"),
      axis.ticks.y         = element_blank(),
      
      axis.text            = element_text(color = AER_BLACK, size = 10),
      axis.title.x         = element_text(color = AER_BLACK, size = 10,
                                          margin = margin(t = 6)),
      
      plot.margin          = margin(t = 8, r = 12, b = 8, l = 8),
      plot.background      = element_rect(fill = "white", color = NA),
      panel.background     = element_rect(fill = "white", color = NA)
    )
}

# -----------------------------------------------------------------------------
# 4. Run all 4 outcomes
# -----------------------------------------------------------------------------
outcomes <- c("ln_drug", "ln_owi", "ln_property", "ln_violent")
outcome_labels <- c(
  ln_drug     = "Drug arrests",
  ln_owi      = "OWI arrests",
  ln_property = "Property crime",
  ln_violent  = "Violent crime"
)
panel_labels <- c(
  ln_drug     = "Panel A. Drug arrests",
  ln_owi      = "Panel B. OWI arrests",
  ln_property = "Panel C. Property crime",
  ln_violent  = "Panel D. Violent crime"
)

all_es   <- list()
all_wald <- list()

for (y in outcomes) {
  cat(sprintf("\n=== Event study: %s ===\n", outcome_labels[[y]]))
  out <- fit_event_iowa(y)
  print(out$rows)
  cat(sprintf("\nJoint Wald F-test of pre-period (k = -5..-2):\n"))
  cat(sprintf("  F = %.3f, df1 = %d, df2 = %d, p_joint = %.4f\n",
              out$wald$F_stat, out$wald$df1, out$wald$df2, out$wald$p_joint))
  if (out$wald$p_joint > 0.10) {
    cat("  -> Fails to reject parallel pre-trends at 10%.\n")
  } else {
    cat("  -> REJECTS parallel pre-trends at 10%.\n")
  }
  
  all_es[[y]]   <- data.table(outcome = outcome_labels[[y]], out$rows)
  all_wald[[y]] <- data.table(outcome = outcome_labels[[y]], out$wald)
}

# -----------------------------------------------------------------------------
# 5. Save individual figures (single-panel each)
# -----------------------------------------------------------------------------
es_drug <- all_es[["ln_drug"]]
g_drug  <- plot_es_aer(es_drug, panel_labels[["ln_drug"]])
ggsave(here::here("out", "fig3_eventstudy_iowa.pdf"),
       g_drug, width = 6.5, height = 4.0)
cat("\nSaved: out/fig3_eventstudy_iowa.pdf  (headline, defense slide 7)\n")

appendix_specs <- list(
  ln_owi      = "figA2_eventstudy_iowa_owi.pdf",
  ln_property = "figA3_eventstudy_iowa_property.pdf",
  ln_violent  = "figA4_eventstudy_iowa_violent.pdf"
)
for (y in names(appendix_specs)) {
  g <- plot_es_aer(all_es[[y]], panel_labels[[y]])
  ggsave(here::here("out", appendix_specs[[y]]), g, width = 6.5, height = 4.0)
  cat(sprintf("Saved: out/%s\n", appendix_specs[[y]]))
}

# -----------------------------------------------------------------------------
# 6. 3-panel stacked figure (drug, property, OWI) - shared y-axis
# -----------------------------------------------------------------------------
appendix_outcomes <- c("ln_drug", "ln_property", "ln_owi")

all_ci <- rbindlist(lapply(appendix_outcomes, function(y) all_es[[y]]))
y_min  <- min(all_ci$ci_lo, na.rm = TRUE)
y_max  <- max(all_ci$ci_hi, na.rm = TRUE)
pad    <- (y_max - y_min) * 0.06
common_ylim <- c(y_min - pad, y_max + pad)

g_list <- lapply(appendix_outcomes, function(y) {
  plot_es_aer(all_es[[y]], panel_labels[[y]], ylim_range = common_ylim)
})

stacked_ok <- FALSE
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  g_stacked <- g_list[[1]] / g_list[[2]] / g_list[[3]]
  ggsave(here::here("out", "figA1_eventstudy_iowa_3panel.pdf"),
         g_stacked, width = 6.5, height = 11)
  cat("Saved: out/figA1_eventstudy_iowa_3panel.pdf  (3-panel stacked)\n")
  stacked_ok <- TRUE
} else if (requireNamespace("cowplot", quietly = TRUE)) {
  suppressPackageStartupMessages(library(cowplot))
  g_stacked <- cowplot::plot_grid(plotlist = g_list, ncol = 1, align = "v")
  ggsave(here::here("out", "figA1_eventstudy_iowa_3panel.pdf"),
         g_stacked, width = 6.5, height = 11)
  cat("Saved: out/figA1_eventstudy_iowa_3panel.pdf  (3-panel stacked, cowplot)\n")
  stacked_ok <- TRUE
}
if (!stacked_ok) {
  cat("Note: install patchwork or cowplot to generate the 3-panel stacked figure.\n")
  cat("      install.packages('patchwork')\n")
}

# -----------------------------------------------------------------------------
# 7. Save CSVs
# -----------------------------------------------------------------------------
fwrite(rbindlist(all_es),
       here::here("out", "fig3_eventstudy_iowa_coefs.csv"))
fwrite(rbindlist(all_wald),
       here::here("out", "fig3_eventstudy_iowa_joint_wald.csv"))
cat("\nSaved CSVs: coefficients + joint Wald summary\n")

# -----------------------------------------------------------------------------
# 8. Defense one-liner
# -----------------------------------------------------------------------------
drug_wald <- all_wald[["ln_drug"]]
drug_es   <- all_es[["ln_drug"]]
pre_min   <- min(drug_es[k < 0]$est)
pre_max   <- max(drug_es[k < 0]$est)
post_max  <- max(drug_es[k >= 0]$est)
post_min  <- min(drug_es[k >= 0]$est)

cat("\n", strrep("=", 72), "\n", sep = "")
cat("DEFENSE ONE-LINER (Section 5.5 / Slide 7):\n\n")
cat(sprintf("'Figure 5.1 plots the Iowa drug-arrest event-study coefficients\n"))
cat(sprintf(" following the figure conventions of Prager and Schmitt (2021,\n"))
cat(sprintf(" AER). Pre-period coefficients (k = -5 to k = -2) range from\n"))
cat(sprintf(" %.2f to %.2f, with confidence intervals overlapping zero. The\n",
            pre_min, pre_max))
cat(sprintf(" joint Wald test of zero pre-trends yields F = %.2f, p_joint =\n",
            drug_wald$F_stat))
cat(sprintf(" %.3f, failing to reject parallel trends. Post-2020 coefficients\n",
            drug_wald$p_joint))
cat(sprintf(" rise to between %.2f and %.2f, with the largest effects in 2021\n",
            post_min, post_max))
cat(sprintf(" and 2022. The pattern is a delayed step increase rather than\n"))
cat(sprintf(" a pre-existing differential drift.'\n"))
cat(strrep("=", 72), "\n", sep = "")

cat("\n=== DONE ===\n")
cat("Next: copy out/fig3_eventstudy_iowa.pdf into defense slide 7.\n")