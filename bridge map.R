# =============================================================================
# bridge_map_publication.R
#
# Publication-quality version of the Iowa-Illinois bridge connectivity map.
# Uses real TIGER 2020 county polygons via {tigris} and {sf}, projected to
# Albers Equal Area Conic (EPSG:5070).
#
# Coexists with bridge_map.R (the schematic fallback). Do not run this if
# the network is restricted or {tigris}/{sf} are unavailable — it will exit
# cleanly with a message; render the schematic version instead.
#
# Inputs (in working directory):
#   - iowa_treated_counties.csv   (8 treated counties + bridge counts)
#   - iowa_il_bridges.csv         (11 bridges with IA/IL endpoint coords)
#
# Outputs:
#   - figures/map_iowa_bridges_publication.pdf  (cairo_pdf)
#   - figures/map_iowa_bridges_publication.png  (300 DPI)
#
# Tested target: R 4.3+, ggplot2 3.5, sf 1.0, tigris 2.1, ggrepel 0.9.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Package availability check (graceful exit if missing)
# -----------------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "ggplot2", "ggrepel", "sf", "tigris")
have <- rownames(installed.packages())
missing_pkgs <- setdiff(required_pkgs, have)
if (length(missing_pkgs) > 0) {
  message("[bridge_map_publication] Missing R packages: ",
          paste(missing_pkgs, collapse = ", "))
  message("  Install with:")
  message("    install.packages(c(\"",
          paste(missing_pkgs, collapse = "\", \""), "\"))")
  message("  Then re-run. The schematic version (bridge_map.R) is the fallback.")
  quit(save = "no", status = 2)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(sf)
  library(tigris)
})

options(tigris_use_cache = TRUE)
sf::sf_use_s2(FALSE)  # planar ops on lon/lat for simple boundary intersection

# -----------------------------------------------------------------------------
# 1. Load CSV inputs
# -----------------------------------------------------------------------------
treated <- read_csv("iowa_treated_counties.csv", show_col_types = FALSE) |>
  mutate(fips = sprintf("%05d", as.integer(fips)))
bridges <- read_csv("iowa_il_bridges.csv",       show_col_types = FALSE)

stopifnot(nrow(treated) == 8L,
          "19115" %in% treated$fips,
          treated$n_bridges[treated$fips == "19115"] == 0,
          treated$n_bridges[treated$fips == "19163"] == 3)

# -----------------------------------------------------------------------------
# 2. Fetch TIGER 2020 polygons + state boundary line
# -----------------------------------------------------------------------------
fetch_or_quit <- function(label, expr) {
  tryCatch(expr,
           error = function(e) {
             message("[bridge_map_publication] TIGER fetch failed for ", label, ":")
             message("  ", conditionMessage(e))
             message("  Most likely cause: no network access to www2.census.gov.")
             message("  Render the schematic instead: source(\"bridge_map.R\")")
             quit(save = "no", status = 3)
           })
}

cat("Fetching TIGER 2020 county polygons for IA + IL ...\n")
ia_co <- fetch_or_quit("IA counties",
                       tigris::counties(state = "IA", year = 2020, cb = TRUE, progress_bar = FALSE))
il_co <- fetch_or_quit("IL counties",
                       tigris::counties(state = "IL", year = 2020, cb = TRUE, progress_bar = FALSE))
us_st <- fetch_or_quit("US states",
                       tigris::states(year = 2020, cb = TRUE, progress_bar = FALSE))

ia_co <- sf::st_as_sf(ia_co)
il_co <- sf::st_as_sf(il_co)
us_st <- sf::st_as_sf(us_st)

ia_poly <- us_st[us_st$STUSPS == "IA", ]
il_poly <- us_st[us_st$STUSPS == "IL", ]

# IA-IL shared boundary (linestring) — proxy for the Mississippi River edge
border_line <- tryCatch(
  sf::st_intersection(sf::st_boundary(ia_poly), sf::st_boundary(il_poly)),
  error = function(e) {
    message("[bridge_map_publication] Border-line intersection failed: ",
            conditionMessage(e))
    message("  Falling back to NULL — the river edge will not be drawn.")
    NULL
  }
)

# -----------------------------------------------------------------------------
# 3. Project everything to Albers Equal Area Conic (EPSG:5070)
# -----------------------------------------------------------------------------
crs_albers <- 5070L
ia_co       <- sf::st_transform(ia_co,       crs_albers)
il_co       <- sf::st_transform(il_co,       crs_albers)
if (!is.null(border_line)) {
  border_line <- sf::st_transform(border_line, crs_albers)
}

# Project bridge endpoints and treated centroids by hand (they're plain data,
# not sf), so geom_segment / geom_point / geom_text_repel can use them.
project_xy <- function(df, lon, lat, suffix = "") {
  pts <- sf::st_as_sf(df[, c(lon, lat)], coords = c(lon, lat), crs = 4326)
  pts <- sf::st_transform(pts, crs_albers)
  xy  <- sf::st_coordinates(pts)
  setNames(as.data.frame(xy), paste0(c("x", "y"), suffix))
}

bridges <- dplyr::bind_cols(
  bridges,
  project_xy(bridges, "ia_lon", "ia_lat", suffix = "_ia"),
  project_xy(bridges, "il_lon", "il_lat", suffix = "_il")
)

treated <- dplyr::bind_cols(
  treated,
  project_xy(treated, "lon", "lat")
)

louisa <- treated |> dplyr::filter(n_bridges == 0)

# -----------------------------------------------------------------------------
# 4. Treatment-intensity factor on the IA county sf
# -----------------------------------------------------------------------------
ia_co <- ia_co |>
  dplyr::mutate(GEOID = as.character(GEOID)) |>
  dplyr::left_join(treated |> dplyr::select(fips, n_bridges),
                   by = c("GEOID" = "fips")) |>
  dplyr::mutate(
    treatment_status = dplyr::case_when(
      !is.na(n_bridges) & n_bridges == 0 ~ "0 — reclassify",
      !is.na(n_bridges) & n_bridges == 1 ~ "1 bridge",
      !is.na(n_bridges) & n_bridges == 2 ~ "2 bridges",
      !is.na(n_bridges) & n_bridges >= 3 ~ "3+ bridges",
      TRUE                               ~ "Control"
    ),
    treatment_status = factor(
      treatment_status,
      levels = c("3+ bridges", "2 bridges", "1 bridge",
                 "0 — reclassify", "Control")
    )
  )

n_treated_matched <- sum(!is.na(ia_co$n_bridges))
stopifnot(n_treated_matched == 8L)

palette_fill <- c(
  "3+ bridges"         = "#6A0F0F",
  "2 bridges"          = "#B91C1C",
  "1 bridge"           = "#E89589",
  "0 — reclassify"= "#888888",
  "Control"            = "#FAF6EC"
)

# -----------------------------------------------------------------------------
# 5. Plot extent: project a lon/lat bbox to Albers and use for coord_sf
# -----------------------------------------------------------------------------
bbox_ll <- sf::st_as_sfc(sf::st_bbox(
  c(xmin = -96.7, xmax = -88.0, ymin = 40.2, ymax = 43.6), crs = 4326))
bbox_albers <- sf::st_bbox(sf::st_transform(bbox_ll, crs_albers))

# State labels: pick lon/lat, project to Albers, take xy
state_label_pts <- project_xy(
  data.frame(lon = c(-93.5, -89.0), lat = c(42.95, 42.95)),
  lon = "lon", lat = "lat"
)

# -----------------------------------------------------------------------------
# 6. Build the plot
# -----------------------------------------------------------------------------
p <- ggplot() +
  # Illinois context fill
  geom_sf(data = il_co, fill = "#EDF2FA", color = "#D8DDE6",
          linewidth = 0.1) +
  
  # Iowa counties filled by treatment intensity
  geom_sf(data = ia_co, aes(fill = treatment_status),
          color = "#BBBBBB", linewidth = 0.25) +
  
  # Iowa-Illinois state boundary (proxy for the Mississippi River edge)
  { if (!is.null(border_line))
    geom_sf(data = border_line, color = "#3E7BAF", linewidth = 1.0,
            alpha = 0.9, inherit.aes = FALSE)
    else NULL } +
  
  # Bridges as red segments
  geom_segment(data = bridges,
               aes(x = x_ia, y = y_ia, xend = x_il, yend = y_il),
               color = "#B91C1C", linewidth = 0.7, alpha = 0.9,
               inherit.aes = FALSE) +
  
  # Louisa: black X marker overlay (no bridge)
  geom_point(data = louisa, aes(x = x, y = y),
             shape = 4, color = "black", size = 5, stroke = 1.6,
             inherit.aes = FALSE) +
  
  # Labels for the 8 treated counties
  geom_text_repel(data = treated,
                  aes(x = x, y = y,
                      label = paste0(name, " (", fips, ")")),
                  size = 3.2, color = "#1A1A1A", family = "serif",
                  segment.color = "#777777", segment.size = 0.3,
                  box.padding = 0.6, point.padding = 0.4,
                  min.segment.length = 0.2, max.overlaps = Inf,
                  seed = 42, inherit.aes = FALSE) +
  
  # State annotations
  annotate("text",
           x = state_label_pts$x[1], y = state_label_pts$y[1],
           label = "IOWA — prohibition state · observed in DiD",
           size = 4.6, fontface = "bold", color = "#333333",
           alpha = 0.7, family = "serif") +
  annotate("text",
           x = state_label_pts$x[2], y = state_label_pts$y[2],
           label = "ILLINOIS — legalized cannabis Jan 1, 2020",
           size = 4.6, fontface = "bold", color = "#333333",
           alpha = 0.7, family = "serif") +
  
  # Scales — only show treated levels in legend; "Control" stays a flat color
  scale_fill_manual(
    values = palette_fill,
    name   = "Bridges to Illinois",
    breaks = c("3+ bridges", "2 bridges", "1 bridge", "0 — reclassify"),
    na.value = "#FAF6EC"
  ) +
  
  coord_sf(crs = crs_albers,
           xlim = c(bbox_albers["xmin"], bbox_albers["xmax"]),
           ylim = c(bbox_albers["ymin"], bbox_albers["ymax"]),
           expand = FALSE) +
  
  labs(
    title = "Treatment intensity by bridge connectivity",
    subtitle = paste0("8 Iowa Mississippi River counties — refined ",
                      "treatment group per Sergey's bridge-connectivity ",
                      "critique. Louisa County has zero river bridge."),
    caption = paste0("County polygons + IA-IL boundary: US Census TIGER 2020. ",
                     "Bridges: Wikipedia 'List of crossings of the Upper ",
                     "Mississippi River' and Iowa DOT records. ",
                     "Projection: Albers Equal Area (EPSG:5070).")
  ) +
  theme_void(base_size = 11, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 14, color = "#1A1A1A",
                                 hjust = 0, margin = margin(b = 4)),
    plot.subtitle = element_text(color = "#444444", size = 10, hjust = 0,
                                 margin = margin(b = 12)),
    plot.caption  = element_text(color = "#666666", size = 7,  hjust = 0,
                                 margin = margin(t = 8)),
    legend.position    = c(0.13, 0.20),
    legend.background  = element_rect(fill = "white", color = "#CCCCCC",
                                      linewidth = 0.3),
    legend.title       = element_text(face = "bold", size = 9,  family = "serif"),
    legend.text        = element_text(size = 8.5, family = "serif"),
    legend.key.height  = unit(6, "mm"),
    plot.margin        = margin(15, 15, 15, 15)
  )

# -----------------------------------------------------------------------------
# 7. Save outputs
# -----------------------------------------------------------------------------
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
pdf_path <- "figures/map_iowa_bridges_publication.pdf"
png_path <- "figures/map_iowa_bridges_publication.png"

ggsave(pdf_path, p, width = 12, height = 8, device = cairo_pdf)
ggsave(png_path, p, width = 12, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 8. One-line summary
# -----------------------------------------------------------------------------
get_n <- function(level) sum(ia_co$treatment_status == level, na.rm = TRUE)
pdf_kb <- round(file.info(pdf_path)$size / 1024, 1)
png_kb <- round(file.info(png_path)$size / 1024, 1)
cat(sprintf(
  "[bridge_map_publication] 12x8in @ 300dpi | PDF %s KB, PNG %s KB | 3+:%d (Scott) 2:%d (Clinton,Lee) 1:%d (Dubuque,Jackson,Muscatine,DesMoines) 0:%d (Louisa—X) Control:%d\n",
  pdf_kb, png_kb,
  get_n("3+ bridges"), get_n("2 bridges"),
  get_n("1 bridge"),   get_n("0 — reclassify"),
  get_n("Control")
))