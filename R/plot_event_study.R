#!/usr/bin/env Rscript
# =============================================================================
# plot_event_study.R — Production-quality event study figure
#
# Run AFTER replication_did.R, which saves event_study_coefs.rds
#
# OUTPUTS:
#   figure_event_study.pdf   (vector, for thesis docx and LaTeX deck)
#   figure_event_study.png   (raster 300 DPI, for slides)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

es <- readRDS("event_study_coefs.rds")

# Use thesis numbers as fallback if you haven't run replication yet
if (!exists("es") || nrow(es) == 0) {
  es <- data.table(
    k = c(-5,-4,-3,-2,-1,0,1,2),
    k_label = c("2015","2016","2017","2018","2019","2020","2021","2022"),
    beta = c(-0.103, 0.018, 0.159, 0.173, 0, 0.089, 0.483, 0.276),
    se   = c(0.153, 0.174, 0.155, 0.205, NA, 0.150, 0.203, 0.176),
    p    = c(0.501, 0.916, 0.308, 0.400, NA, 0.555, 0.018, 0.116)
  )
  es[, ci_lo := beta - 1.96 * se]
  es[, ci_hi := beta + 1.96 * se]
}

# Static DiD value for reference line
ATT_STATIC <- 0.233
ATT_SE <- 0.103
ATT_P <- 0.024

# Pre-trend Wald
WALD_F <- 0.86
WALD_P <- 0.486

# Year labels with k below in subscript
es[, x_label := paste0(k_label, "\n(k=", k, ")")]
es[, is_ref := k == -1]
es[, is_post := k >= 0]

# Color palette — neutral and print-friendly
col_pre <- "#2c3e50"      # dark slate
col_post <- "#1f4e79"     # dark blue
col_ref <- "#888888"      # gray
col_att <- "#c0392b"      # dark red
col_att_band <- "#fadbd8" # light red

p <- ggplot(es, aes(x = k, y = beta)) +
  
  # Background ATT band (between pre/post split)
  annotate("rect",
           xmin = -0.5, xmax = 2.5,
           ymin = ATT_STATIC - 1.96 * ATT_SE,
           ymax = ATT_STATIC + 1.96 * ATT_SE,
           fill = col_att_band, alpha = 0.4) +
  
  # ATT reference horizontal line
  geom_hline(yintercept = ATT_STATIC, color = col_att,
             linetype = "solid", linewidth = 0.6, alpha = 0.7) +
  
  # Zero reference line
  geom_hline(yintercept = 0, color = "black",
             linetype = "solid", linewidth = 0.4) +
  
  # Vertical line at treatment (between k=-1 and k=0)
  geom_vline(xintercept = -0.5, color = "#888888",
             linetype = "dashed", linewidth = 0.5) +
  
  # CI error bars
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi,
                    color = ifelse(is_post, "post", "pre")),
                width = 0.15, linewidth = 0.7, na.rm = TRUE) +
  
  # Coefficient points
  geom_point(aes(color = ifelse(is_ref, "ref",
                                ifelse(is_post, "post", "pre")),
                 shape = is_ref),
             size = 3, stroke = 1.2) +
  
  # Reference point (k=-1) styling
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16), guide = "none") +
  
  # Color manual
  scale_color_manual(values = c("pre" = col_pre, "post" = col_post, "ref" = col_ref),
                     guide = "none") +
  
  # Annotations
  annotate("text", x = -0.4, y = max(es$ci_hi, na.rm=TRUE) * 0.95,
           label = "Jan 2020\nIL legalization",
           hjust = 0, vjust = 1, size = 3, color = "#666666",
           fontface = "italic") +
  
  annotate("text", x = 1.7, y = ATT_STATIC + 0.04,
           label = paste0("Static ATT = ", sprintf("%.3f", ATT_STATIC)),
           color = col_att, size = 3.2, hjust = 0.5,
           fontface = "bold") +
  
  # Pre-trend test box (upper left)
  annotate("rect", xmin = -5.4, xmax = -1.6,
           ymin = max(es$ci_hi, na.rm = TRUE) - 0.08,
           ymax = max(es$ci_hi, na.rm = TRUE) + 0.02,
           fill = "#f5f5f5", color = "#bbbbbb", linewidth = 0.3) +
  annotate("text", x = -5.3, y = max(es$ci_hi, na.rm = TRUE) - 0.01,
           label = paste0("Pre-trend joint Wald: F = ",
                          sprintf("%.2f", WALD_F),
                          ", p = ", sprintf("%.3f", WALD_P)),
           hjust = 0, vjust = 1, size = 3, family = "sans") +
  
  # Post-2020 ATT box (lower right)
  annotate("rect", xmin = 0.6, xmax = 2.4,
           ymin = min(es$ci_lo, na.rm = TRUE) - 0.02,
           ymax = min(es$ci_lo, na.rm = TRUE) + 0.10,
           fill = "#f5f5f5", color = "#bbbbbb", linewidth = 0.3) +
  annotate("text", x = 0.7, y = min(es$ci_lo, na.rm = TRUE) + 0.085,
           label = paste0("Post-2020 ATT = ", sprintf("%.3f", ATT_STATIC)),
           hjust = 0, vjust = 1, size = 3, fontface = "bold") +
  annotate("text", x = 0.7, y = min(es$ci_lo, na.rm = TRUE) + 0.045,
           label = paste0("HC1 SE = ", sprintf("%.3f", ATT_SE),
                          ", p = ", sprintf("%.3f", ATT_P)),
           hjust = 0, vjust = 1, size = 2.8, color = "#444444") +
  annotate("text", x = 0.7, y = min(es$ci_lo, na.rm = TRUE) + 0.005,
           label = paste0("95% CI: [",
                          sprintf("%.2f", ATT_STATIC - 1.96 * ATT_SE), ", ",
                          sprintf("%.2f", ATT_STATIC + 1.96 * ATT_SE), "]"),
           hjust = 0, vjust = 1, size = 2.8, color = "#444444") +
  
  # Scales
  scale_x_continuous(breaks = -5:2,
                     labels = paste0(es$k_label, "\n(", es$k, ")")) +
  scale_y_continuous(breaks = seq(-0.6, 1.0, 0.2)) +
  
  # Labels
  labs(
    title = "Event study: drug-enforcement arrests in Iowa bridge counties",
    subtitle = "Relative to interior counties, 2015-2022",
    x = "Years relative to January 2020",
    y = expression("Event-study coefficient " * beta[k] * " on log(drug arrests + 1)"),
    caption = paste0(
      "Notes: 7-bridge specification on the predefined Iowa DOT bridge inventory. ",
      "Treated: Clinton, Des Moines, Dubuque, Jackson, Lee, Muscatine, Scott. ",
      "Control: 92 Iowa interior counties. ",
      "Cluster-robust SE at county level (HC1). Reference period k = -1 (2019). ",
      "Pre-period joint Wald test F = ", sprintf("%.2f", WALD_F),
      ", p = ", sprintf("%.3f", WALD_P), " (fails to reject parallel trends)."
    )
  ) +
  
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "#555555", hjust = 0,
                                 margin = margin(b = 8)),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    axis.text = element_text(size = 9, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.3),
    plot.caption = element_text(size = 7.5, color = "#666666",
                                hjust = 0, margin = margin(t = 12),
                                lineheight = 1.15),
    plot.margin = margin(15, 20, 12, 15)
  )

# Save vector PDF for thesis
ggsave("figure_event_study.pdf", p, width = 10, height = 6.2,
       units = "in", device = cairo_pdf)

# Save PNG for slides (300 DPI)
ggsave("figure_event_study.png", p, width = 10, height = 6.2,
       units = "in", dpi = 300)

cat("Saved figure_event_study.pdf (vector, for thesis)\n")
cat("Saved figure_event_study.png (300 DPI, for slides)\n")
