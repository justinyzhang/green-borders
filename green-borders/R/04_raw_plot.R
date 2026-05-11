# =============================================================================
# 04_raw_plot.R   (4-state version, drug headline)
#
# Figure 2: raw mean drug arrest rate by IL-border status, IA + IN + WI,
# 2015-2022. Headline figure for Chapter 5.
#
# Outputs:
#   out/fig2_raw_means.pdf       (drug headline)
#   out/figA1_raw_means_all.pdf  (4-panel appendix)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
})

panel <- readRDS(here::here("data", "processed", "panel_county_year.rds"))
setDT(panel)
panel_main <- panel[treatment %in% c("IL_border", "interior")]

outcomes <- c("drug", "owi", "property", "violent")
outcome_labels <- c(
  drug     = "Drug arrests",
  owi      = "OWI arrests",
  property = "Property crime",
  violent  = "Violent crime"
)

agg_long <- rbindlist(lapply(outcomes, function(o) {
  rate_col <- paste0(o, "_rate")
  panel_main[, .(year, outcome = outcome_labels[[o]],
                 rate = get(rate_col), is_il_border)]
}))[, .(mean_rate = mean(rate, na.rm = TRUE)),
    by = .(year, outcome,
           group = ifelse(is_il_border == 1,
                          "IL-border (n = 27)", "Interior (n = 236)"))]

agg_drug <- agg_long[outcome == "Drug arrests"]

g_drug <- ggplot(agg_drug, aes(year, mean_rate,
                               color = group, shape = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 2019.5, linetype = "dashed", color = "gray30") +
  annotate("text", x = 2019.6,
           y = max(agg_drug$mean_rate, na.rm = TRUE) * 0.98,
           label = "IL legalization", hjust = 0, size = 3.3, color = "gray20") +
  scale_color_manual(values = c("IL-border (n = 27)" = "#1f78b4",
                                "Interior (n = 236)" = "#7a7a7a")) +
  scale_shape_manual(values = c("IL-border (n = 27)" = 16,
                                "Interior (n = 236)" = 17)) +
  scale_x_continuous(breaks = 2015:2022) +
  labs(x = NULL, y = "Drug arrests per 100,000",
       color = NULL, shape = NULL,
       title = "Drug arrest rate, IL-border vs interior counties in IA + IN + WI",
       subtitle = "Annual county means; vertical line = January 2020 IL legalization") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

dir.create(here::here("out"), showWarnings = FALSE, recursive = TRUE)
ggsave(here::here("out", "fig2_raw_means.pdf"), g_drug, width = 7.5, height = 4.5)
cat("Wrote out/fig2_raw_means.pdf\n")

agg_long[, outcome := factor(outcome, levels = c("Drug arrests", "OWI arrests",
                                                 "Property crime", "Violent crime"))]
g_all <- ggplot(agg_long, aes(year, mean_rate,
                              color = group, shape = group, group = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = 2019.5, linetype = "dashed", color = "gray30") +
  facet_wrap(~ outcome, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("IL-border (n = 27)" = "#1f78b4",
                                "Interior (n = 236)" = "#7a7a7a")) +
  scale_shape_manual(values = c("IL-border (n = 27)" = 16,
                                "Interior (n = 236)" = 17)) +
  scale_x_continuous(breaks = c(2015, 2018, 2021)) +
  labs(x = NULL, y = "Arrests per 100,000",
       color = NULL, shape = NULL,
       title = "Crime outcomes by IL-border status, IA + IN + WI, 2015-2022") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave(here::here("out", "figA1_raw_means_all.pdf"), g_all, width = 8, height = 6)
cat("Wrote out/figA1_raw_means_all.pdf\n")
