# ============================================================================
# R/10_permutation_test.R
# ----------------------------------------------------------------------------
# Permutation / placebo test for the Iowa DiD headline result.
#
# Logic: Randomly reassign "treatment" to 8 of the 99 Iowa counties,
#        re-estimate the baseline DiD, and repeat B = 1000 times.
#        The real β = 0.422 should lie in the extreme tail of the
#        placebo distribution if the result is not driven by chance.
#
# This is the most stringent inference test available because it does NOT
# rely on cluster-SE asymptotics at all. Particularly important given only
# 8 treated clusters (Cameron & Miller 2015 small-N concern).
#
# Defense framing: addresses Concern #1 in the playbook
#   "8 cluster cluster-robust SE 不可靠" → "I ran permutation test, real β
#    is in top X% of 1000 placebo assignments, p_perm = X.XX"
#
# Inputs:  data/processed/panel_county_year.rds
# Outputs: out/permutation_test_results.rds  (placebo distribution + p)
#          out/permutation_test_plot.pdf      (histogram visualization)
#          out/permutation_test_summary.txt   (defense one-liner)
#
# Wall-clock: ~2-5 min on a 2020-era laptop for B = 1000 iterations
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("data.table", "fixest", "ggplot2")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(data.table)
library(fixest)
library(ggplot2)

dir.create("out", showWarnings = FALSE, recursive = TRUE)

set.seed(42)
B <- 1000                       # number of permutations
N_TREATED <- 8                  # match the true treated cluster count

# ---- 1. Load panel ---------------------------------------------------------

panel <- readRDS("data/processed/panel_county_year.rds")
panel <- as.data.table(panel)

state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
ia_panel <- panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia")]
stopifnot(nrow(ia_panel) > 0)

cat(sprintf("Iowa panel: %d county-year obs\n", nrow(ia_panel)))
all_ia_counties <- unique(ia_panel$county_fips)
cat(sprintf("Total Iowa counties: %d\n", length(all_ia_counties)))

# Verify real treated set has 8 counties
real_treated <- unique(ia_panel[border == 1, county_fips])
cat(sprintf("Real treated (border) counties: %d\n", length(real_treated)))
stopifnot(length(real_treated) == N_TREATED)

# ---- 2. Get the REAL β as benchmark ---------------------------------------

cat("\n=== Real (true) Iowa DiD coefficient ===\n")

if (!"border_x_post" %in% names(ia_panel)) {
  ia_panel[, border_x_post := border * post]
}

m_real <- feols(
  ln_drug ~ border_x_post + ln_pop + ln_inc | county_fips + year,
  data = ia_panel,
  cluster = ~county_fips
)

real_beta <- coef(m_real)["border_x_post"]
real_p_hc1 <- pvalue(m_real)["border_x_post"]
cat(sprintf("Real β = %.4f (p_HC1 = %.4f)\n", real_beta, real_p_hc1))

# Sanity check: should be ~0.42
if (abs(real_beta - 0.42) > 0.05) {
  warning(sprintf("Real β = %.3f differs from expected 0.42 — check panel build", real_beta))
}

# ---- 3. Permutation loop ---------------------------------------------------

cat(sprintf("\n=== Running %d placebo permutations ===\n", B))
cat("Each iteration: randomly assign 8 of 99 IA counties as 'treated',\n")
cat("then re-estimate the DiD coefficient on log(drug arrests + 1).\n\n")

placebo_betas <- numeric(B)
progress_step <- B / 20  # print every 5%

for (b in seq_len(B)) {
  # Random reassignment of 8 fake treated counties
  fake_treated <- sample(all_ia_counties, N_TREATED)
  
  ia_panel[, fake_border := as.integer(county_fips %in% fake_treated)]
  ia_panel[, fake_border_x_post := fake_border * post]
  
  # Fit the same spec with fake treatment
  m_fake <- tryCatch(
    feols(
      ln_drug ~ fake_border_x_post + ln_pop + ln_inc | county_fips + year,
      data = ia_panel,
      cluster = ~county_fips,
      warn = FALSE, notes = FALSE
    ),
    error = function(e) NULL
  )
  
  placebo_betas[b] <- if (!is.null(m_fake)) coef(m_fake)["fake_border_x_post"] else NA_real_
  
  if (b %% progress_step == 0) {
    cat(sprintf("  ... %d / %d (%.0f%%)\n", b, B, 100 * b / B))
  }
}

# Drop any failed iterations
placebo_betas <- placebo_betas[!is.na(placebo_betas)]
B_actual <- length(placebo_betas)

# ---- 4. Compute permutation p-values ---------------------------------------

# Two-sided p: fraction of placebo |β| >= |real β|
p_perm_2sided <- mean(abs(placebo_betas) >= abs(real_beta))
# One-sided p: fraction of placebo β >= real β
p_perm_1sided <- mean(placebo_betas >= real_beta)

# Where does real_beta sit in the placebo distribution?
percentile_real <- mean(placebo_betas < real_beta) * 100

cat("\n=== Results ===\n")
cat(sprintf("Real β:                  %+.4f\n", real_beta))
cat(sprintf("Placebo β mean:          %+.4f\n", mean(placebo_betas)))
cat(sprintf("Placebo β median:        %+.4f\n", median(placebo_betas)))
cat(sprintf("Placebo β SD:            %.4f\n", sd(placebo_betas)))
cat(sprintf("Placebo β 2.5%% / 97.5%%:  [%.4f, %.4f]\n",
            quantile(placebo_betas, 0.025), quantile(placebo_betas, 0.975)))
cat(sprintf("\nReal β percentile in placebo distribution: %.1f%%\n", percentile_real))
cat(sprintf("Permutation p (2-sided): %.4f\n", p_perm_2sided))
cat(sprintf("Permutation p (1-sided): %.4f\n", p_perm_1sided))
cat(sprintf("HC1 p (for comparison):  %.4f\n", real_p_hc1))

# ---- 5. Save outputs -------------------------------------------------------

results <- list(
  real_beta = real_beta,
  real_p_hc1 = real_p_hc1,
  placebo_betas = placebo_betas,
  B_actual = B_actual,
  p_perm_2sided = p_perm_2sided,
  p_perm_1sided = p_perm_1sided,
  percentile_real = percentile_real
)
saveRDS(results, "out/permutation_test_results.rds")

# Plot the placebo distribution + real β vertical line
p_plot <- ggplot(data.frame(beta = placebo_betas), aes(x = beta)) +
  geom_histogram(bins = 40, fill = "grey75", color = "white") +
  geom_vline(xintercept = real_beta, color = "#B91C1C", linewidth = 1.0) +
  annotate("text", x = real_beta, y = Inf,
           label = sprintf(" Real β = %.3f\n p_perm = %.3f", real_beta, p_perm_2sided),
           hjust = 0, vjust = 2, color = "#B91C1C", size = 4.2) +
  labs(
    title    = "Permutation test: real β vs 1000 placebo assignments",
    subtitle = sprintf("Iowa drug arrests, B = %d random reassignments of 8 fake treated counties",
                       B_actual),
    x = expression(paste("Placebo coefficient ", beta)),
    y = "Frequency",
    caption = "Cameron & Miller (2015) inference robust to small-cluster asymptotics"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"))

ggsave("out/permutation_test_plot.pdf", p_plot, width = 8, height = 5)

# Human-readable summary
sink("out/permutation_test_summary.txt")
cat("Permutation Test — Iowa drug arrest DiD\n")
cat("=========================================\n\n")
cat(sprintf("Number of permutations:        B = %d\n", B_actual))
cat(sprintf("Real β:                        %+.4f\n", real_beta))
cat(sprintf("Real β percentile in placebo:  %.1f%%\n", percentile_real))
cat(sprintf("Permutation p (two-sided):     %.4f\n", p_perm_2sided))
cat(sprintf("Permutation p (one-sided):     %.4f\n", p_perm_1sided))
cat(sprintf("HC1 p (cluster-robust SE):     %.4f\n", real_p_hc1))
cat("\nInterpretation:\n")
cat(" - p_perm tests whether the real β is unusually large compared\n")
cat("   to random reassignments of the treatment label.\n")
cat(" - This inference does NOT rely on cluster-SE asymptotics and is\n")
cat("   recommended by Cameron & Miller (2015) when N_clusters is small.\n")
cat(" - p_perm < 0.05 → real β is in the top 5%% of placebo distribution\n")
cat("   → rejection cannot be attributed to small-N inference artifacts.\n")
sink()

cat("\n=== Saved ===\n")
cat("  out/permutation_test_results.rds  (full placebo distribution)\n")
cat("  out/permutation_test_plot.pdf     (histogram visualization)\n")
cat("  out/permutation_test_summary.txt  (human-readable)\n")

# ---- 6. Defense one-liner --------------------------------------------------

cat("\n",
    strrep("─", 70), "\n",
    "DEFENSE ONE-LINER:\n",
    sprintf("'I ran a permutation test with B = %d random reassignments\n", B_actual),
    sprintf(" of the 8 treated counties to other Iowa counties. The real\n"),
    sprintf(" coefficient β = %.3f lies at the %.0fth percentile of the\n", real_beta, percentile_real),
    sprintf(" placebo distribution, with a two-sided permutation p-value of\n"),
    sprintf(" %.3f. This inference does not depend on cluster-SE asymptotics\n", p_perm_2sided),
    " and directly addresses the small-treated-cluster concern.'\n",
    strrep("─", 70), "\n",
    sep="")

cat("\nDone. Add permutation row to Table 5.2 in thesis Ch5.\n")
