#!/usr/bin/env Rscript
# =============================================================================
# plot_il_heatmap.R — Production-quality Illinois dispensary spatial map
#
# Shows Illinois adult-use dispensary concentration with Iowa 7 bridge
# counties highlighted, demonstrating the supply-side cross-border pattern.
#
# REQUIRED INPUTS:
#   - dispensaries_clean.csv  (dispensary locations with lat/lon)
#
# REQUIRED R PACKAGES:
#   install.packages(c("sf","ggplot2","data.table","tigris","viridis"))
#
# OUTPUTS:
#   figure_il_heatmap.pdf    (vector, for thesis)
#   figure_il_heatmap.png    (300 DPI for slides)
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(data.table)
  library(viridis)
})

# =============================================================================
# LOAD DISPENSARY DATA
# =============================================================================
disp <- fread("dispensaries_clean.csv")
cat("Dispensary records:", nrow(disp), "\n")
cat("Columns:", paste(names(disp), collapse = ", "), "\n")

# Adjust column names if needed (your CSV may have different naming)
if (!all(c("lat", "lon") %in% names(disp))) {
  if ("latitude" %in% names(disp)) setnames(disp, "latitude", "lat")
  if ("longitude" %in% names(disp)) setnames(disp, "longitude", "lon")
}

# Sanity-check coordinates fall in Illinois region
disp <- disp[!is.na(lat) & !is.na(lon)]
disp <- disp[lat > 36 & lat < 43 & lon > -92 & lon < -87]
cat("Dispensaries within Illinois bounding box:", nrow(disp), "\n")

# =============================================================================
# LOAD STATE & COUNTY GEOMETRIES
# =============================================================================
# Option A: Use tigris (downloads from Census, requires internet)
if (requireNamespace("tigris", quietly = TRUE)) {
  library(tigris)
  options(tigris_use_cache = TRUE)
  
  # IL and IA counties
  il_counties <- counties(state = "IL", year = 2020, cb = TRUE, progress_bar = FALSE)
  ia_counties <- counties(state = "IA", year = 2020, cb = TRUE, progress_bar = FALSE)
  states <- states(year = 2020, cb = TRUE, progress_bar = FALSE)
  
  # Filter neighboring states for map context
  states_show <- states[states$STUSPS %in% c("IL", "IA", "WI", "IN", "MO", "KY"), ]
  states_show <- st_transform(states_show, 5070)  # Albers equal-area
  il_counties <- st_transform(il_counties, 5070)
  ia_counties <- st_transform(ia_counties, 5070)
} else {
  stop("tigris package required for downloading geographies. Install with install.packages('tigris')")
}

# =============================================================================
# IDENTIFY 7 IOWA BRIDGE COUNTIES
# =============================================================================
BRIDGE_7 <- c("19045","19057","19061","19097","19111","19139","19163")
ia_counties$is_bridge <- ia_counties$GEOID %in% BRIDGE_7

bridge_geom <- ia_counties[ia_counties$is_bridge, ]
ia_interior <- ia_counties[!ia_counties$is_bridge, ]

# =============================================================================
# AGGREGATE DISPENSARIES TO IL COUNTY
# =============================================================================
disp_sf <- st_as_sf(disp, coords = c("lon", "lat"), crs = 4326)
disp_sf <- st_transform(disp_sf, 5070)

# Spatial join: dispensary → IL county
disp_with_county <- st_join(disp_sf, il_counties[, c("GEOID","NAME")])
disp_count <- as.data.table(st_drop_geometry(disp_with_county))[, .N, by = GEOID]

il_counties_dt <- as.data.table(il_counties)
il_counties_dt <- merge(il_counties_dt, disp_count, by = "GEOID", all.x = TRUE)
il_counties_dt[is.na(N), N := 0]
il_counties$disp_count <- il_counties_dt$N[match(il_counties$GEOID, il_counties_dt$GEOID)]

# Per-100k normalization (need population — fallback to raw count if not available)
if ("ALAND" %in% names(il_counties)) {
  # Use county area as rough density proxy (per square km, not per 100k)
  # If you have population data, replace this with disp_count / pop * 100000
  il_counties$disp_per_km2 <- il_counties$disp_count / (as.numeric(il_counties$ALAND) / 1e6)
}

# =============================================================================
# BUILD PLOT
# =============================================================================
# Bounding box for zoom: IL + eastern IA + bits of WI/IN/MO
bbox <- st_bbox(c(xmin = -91.8, ymin = 36.8, xmax = -86.8, ymax = 42.6),
                crs = st_crs(4326))
bbox_proj <- st_bbox(st_transform(st_as_sfc(bbox), 5070))

p <- ggplot() +
  
  # Background states (light gray)
  geom_sf(data = states_show, fill = "#fafafa", color = "#cccccc", linewidth = 0.3) +
  
  # Illinois counties choropleth (dispensary count)
  geom_sf(data = il_counties, aes(fill = disp_count),
          color = "white", linewidth = 0.15) +
  
  # Iowa interior counties (outline only)
  geom_sf(data = ia_interior, fill = "#f5f5f5", color = "#aaaaaa",
          linewidth = 0.2) +
  
  # Iowa 7 BRIDGE counties (highlighted)
  geom_sf(data = bridge_geom, fill = "#fff4e6", color = "#c0392b",
          linewidth = 0.8) +
  
  # State boundaries (thick)
  geom_sf(data = states_show, fill = NA, color = "black", linewidth = 0.5) +
  
  # Dispensary points
  geom_sf(data = disp_sf, color = "#2c3e50", size = 0.5, alpha = 0.6) +
  
  # Bridge county labels
  geom_sf_text(data = bridge_geom,
               aes(label = stringr::str_to_title(NAME)),
               size = 2.5, color = "#7b241c", fontface = "bold",
               nudge_x = -20000) +
  
  # Color scale
  scale_fill_viridis_c(
    name = "IL dispensary\ncount (per county)",
    option = "plasma",
    direction = -1,
    breaks = c(0, 5, 10, 20),
    na.value = "white"
  ) +
  
  # Zoom to area of interest
  coord_sf(
    xlim = c(bbox_proj["xmin"], bbox_proj["xmax"]),
    ylim = c(bbox_proj["ymin"], bbox_proj["ymax"]),
    expand = FALSE
  ) +
  
  # Labels
  labs(
    title = "Illinois adult-use cannabis dispensaries cluster near prohibition borders",
    subtitle = paste0(nrow(disp), " licensed adult-use dispensaries through 2024 (IDFPR registry)"),
    caption = paste0(
      "Notes: Illinois counties colored by cumulative adult-use dispensary licenses 2019-2024. ",
      "Iowa 7 bridge counties outlined in red: Clinton, Des Moines, Dubuque, Jackson, Lee, ",
      "Muscatine, Scott (Iowa DOT bridge inventory). Iowa interior counties shown in light gray. ",
      "Bordering states (WI, IN, KY, MO) in light background.\n",
      "Sources: IDFPR adult-use dispensary registry; US Census TIGER/Line 2020 county geometry; ",
      "Iowa DOT bridge inventory."
    )
  ) +
  
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#555555", hjust = 0,
                                  margin = margin(b = 10)),
    legend.position = c(0.92, 0.18),
    legend.background = element_rect(fill = "white", color = "#cccccc",
                                      linewidth = 0.3),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.caption = element_text(size = 7.5, color = "#666666",
                                 hjust = 0, margin = margin(t = 8),
                                 lineheight = 1.2),
    plot.margin = margin(15, 15, 12, 15)
  )

# Save
ggsave("figure_il_heatmap.pdf", p, width = 10, height = 7,
       units = "in", device = cairo_pdf)
ggsave("figure_il_heatmap.png", p, width = 10, height = 7,
       units = "in", dpi = 300)

cat("Saved figure_il_heatmap.pdf (vector, for thesis)\n")
cat("Saved figure_il_heatmap.png (300 DPI, for slides)\n")
