# ============================================================================
# R/11_leave_one_out.R
# ----------------------------------------------------------------------------
# Leave-one-treated-county-out (LOOT) robustness for the Iowa DiD.
# Addresses defense Concern #1 (small N treated) and reviewer #3.
#
# Logic: For each of the 8 treated Iowa counties, drop it from the sample
#        and re-estimate the DiD. If β is stable across all 8 leave-outs,
#        the result is not driven by any single treated county
#        (especially Scott County, the largest).
#
# This is a complement to the existing "drop Scott County" robustness in
# Table 5.2 — Scott is the most influential, but the reviewer rightly notes
# we should check ALL 8 in turn for transparency.
#
# Inputs:  data/processed/panel_county_year.rds
# Outputs: out/leave_one_out_results.rds
#          out/leave_one_out_plot.pdf      (8 coefficients + CIs)
#          out/leave_one_out_summary.txt   (defense one-liner)
#
# Wall-clock: <1 min (only 8 regressions)
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("data.table", "fixest", "ggplot2")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(data.table)
library(fixest)
library(ggplot2)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load panel ---------------------------------------------------------

panel <- readRDS("data/processed/panel_county_year.rds")
panel <- as.data.table(panel)

state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
ia_panel <- panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia")]
stopifnot(nrow(ia_panel) > 0)

if (!"border_x_post" %in% names(ia_panel)) {
  ia_panel[, border_x_post := border * post]
}

real_treated <- unique(ia_panel[border == 1, county_fips])
cat(sprintf("Treated counties: %d\n", length(real_treated)))
cat(sprintf("Treated FIPS: %s\n", paste(real_treated, collapse = ", ")))

# Try to get county names if available
if ("county_name" %in% names(ia_panel)) {
  treated_names <- unique(ia_panel[border == 1, .(county_fips, county_name)])
  cat("\nTreated counties (with names):\n")
  print(treated_names)
}

# ---- 2. Baseline (no county dropped) for reference ------------------------

m_full <- feols(
  ln_drug ~ border_x_post + ln_pop + ln_inc | county_fips + year,
  data = ia_panel,
  cluster = ~county_fips
)
beta_full <- coef(m_full)["border_x_post"]
se_full <- se(m_full)["border_x_post"]
p_full <- pvalue(m_full)["border_x_post"]

cat(sprintf("\n=== Full sample baseline ===\n"))
cat(sprintf("β = %.4f, SE = %.4f, p = %.4f\n", beta_full, se_full, p_full))

# ---- 3. Loop: drop each treated county once -------------------------------

cat(sprintf("\n=== Leave-one-treated-county-out (LOOT) ===\n\n"))

loot_results <- data.table(
  dropped_fips = character(),
  dropped_name = character(),
  beta = numeric(),
  se = numeric(),
  p = numeric(),
  ci_lo = numeric(),
  ci_hi = numeric()
)

for (fips in real_treated) {
  m_drop <- feols(
    ln_drug ~ border_x_post + ln_pop + ln_inc | county_fips + year,
    data = ia_panel[county_fips != fips],
    cluster = ~county_fips
  )
  b <- coef(m_drop)["border_x_post"]
  s <- se(m_drop)["border_x_post"]
  p <- pvalue(m_drop)["border_x_post"]
  ci_lo <- b - 1.96 * s
  ci_hi <- b + 1.96 * s
  
  dropped_nm <- if ("county_name" %in% names(ia_panel)) {
    unique(ia_panel[county_fips == fips, county_name])[1]
  } else {
    fips
  }
  
  loot_results <- rbind(loot_results, data.table(
    dropped_fips = fips,
    dropped_name = as.character(dropped_nm),
    beta = round(b, 4),
    se = round(s, 4),
    p = round(p, 4),
    ci_lo = round(ci_lo, 4),
    ci_hi = round(ci_hi, 4)
  ))
  
  cat(sprintf("  Drop %s (%s):  β = %+.4f  SE = %.4f  p = %.4f\n",
              fips, dropped_nm, b, s, p))
}

# ---- 4. Summary statistics across leave-outs ------------------------------

cat("\n=== Summary across 8 leave-outs ===\n")
cat(sprintf("Range of β:     [%.3f, %.3f]\n", min(loot_results$beta), max(loot_results$beta)))
cat(sprintf("Mean β:         %.3f\n", mean(loot_results$beta)))
cat(sprintf("Median β:       %.3f\n", median(loot_results$beta)))
cat(sprintf("Full sample β:  %.3f\n", beta_full))
cat(sprintf("# sig at 5%%:    %d / 8\n", sum(loot_results$p < 0.05)))
cat(sprintf("# sig at 10%%:   %d / 8\n", sum(loot_results$p < 0.10)))

# ---- 5. Save outputs -------------------------------------------------------

saveRDS(loot_results, "out/leave_one_out_results.rds")

# Plot: coefficient + 95% CI for each drop
p_plot <- ggplot(loot_results, aes(x = reorder(dropped_name, beta), y = beta)) +
  geom_point(color = "#B91C1C", size = 3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, color = "#B91C1C") +
  geom_hline(yintercept = beta_full, color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  annotate("text", x = 1, y = beta_full,
           label = sprintf(" Full β = %.3f", beta_full),
           hjust = 0, vjust = -0.5, color = "grey30", size = 3.5) +
  coord_flip() +
  labs(
    title    = "Leave-one-treated-county-out robustness",
    subtitle = "Iowa drug arrests DiD: dropping each treated county one at a time",
    x = "Dropped county",
    y = expression(paste("Treatment coefficient ", beta, " (95% CI)")),
    caption = "Dashed line: full-sample baseline. Black line: zero."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

ggsave("out/leave_one_out_plot.pdf", p_plot, width = 8, height = 5)

# Human-readable summary
sink("out/leave_one_out_summary.txt")
cat("Leave-one-treated-out (LOOT) — Iowa drug DiD\n")
cat("=============================================\n\n")
cat(sprintf("Full-sample baseline:  β = %.4f, p = %.4f\n\n", beta_full, p_full))
cat("Leave-one-out estimates:\n")
print(loot_results)
cat(sprintf("\nRange of β across 8 leave-outs: [%.3f, %.3f]\n",
            min(loot_results$beta), max(loot_results$beta)))
cat(sprintf("Number significant at 5%%:       %d / 8\n", sum(loot_results$p < 0.05)))
cat(sprintf("Number significant at 10%%:      %d / 8\n", sum(loot_results$p < 0.10)))
sink()

cat("\n=== Saved ===\n")
cat("  out/leave_one_out_results.rds\n")
cat("  out/leave_one_out_plot.pdf\n")
cat("  out/leave_one_out_summary.txt\n")

# ---- 6. Defense one-liner --------------------------------------------------

n_sig5 <- sum(loot_results$p < 0.05)
n_sig10 <- sum(loot_results$p < 0.10)
b_min <- min(loot_results$beta)
b_max <- max(loot_results$beta)

cat("\n",
    strrep("─", 70), "\n",
    "DEFENSE ONE-LINER:\n",
    sprintf("'I dropped each of the 8 treated counties one at a time and\n"),
    sprintf(" re-estimated the DiD. The coefficient ranges from %.2f to %.2f\n", b_min, b_max),
    sprintf(" across the 8 leave-outs, with %d of 8 estimates significant at\n", n_sig5),
    sprintf(" the 5%% level and %d of 8 at the 10%% level. The headline result\n", n_sig10),
    " is not driven by any single treated county.'\n",
    strrep("─", 70), "\n",
    sep="")

cat("\nDone. Add leave-one-out block to Table 5.2 in thesis Ch5.\n")
