# ============================================================================
# R/13_county_maps.R
# ----------------------------------------------------------------------------
# Two county-level maps for the Green Borders thesis:
#
#  Map A: Illinois + all 5 prohibition neighbors (WI, IA, MO, IN, KY)
#         - IL counties shaded by dispensary density
#         - State borders highlighted
#         - Border15 counties marked
#
#  Map B: Illinois + Iowa only
#         - IL counties shaded by dispensary density
#         - 8 IA treated counties highlighted in red
#         - 91 IA control counties in light grey
#         - Mississippi River boundary visible
#
# Inputs:
#   - data/processed/il_cross_section.rds (or similar)  -> for IL dispensary data
#     Required columns: county_fips (5-digit), EntryPer100k, Border15
#   - data/processed/panel_county_year.rds  -> for IA treated indicator
#     Required columns: county_fips, state_name, border
#
# Outputs:
#   - out/map_A_all_borders.pdf  +  out/map_A_all_borders.png
#   - out/map_B_iowa_focus.pdf   +  out/map_B_iowa_focus.png
#
# Wall-clock: ~30-60 sec (mostly TIGER download on first run; cached after)
# ============================================================================

# ---- 0. Setup --------------------------------------------------------------

pkgs <- c("data.table", "sf", "tigris", "ggplot2", "viridis", "scales")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(data.table)
library(sf)
library(tigris)
library(ggplot2)
library(viridis)
library(scales)

# Cache tigris shapefiles for repeat runs
options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

dir.create("out", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load county geometry from TIGER -----------------------------------

cat("Pulling TIGER county geometry for IL + 5 neighbors...\n")

# 6 states: IL + 5 prohibition neighbors
state_fips <- c("17", "19", "55", "29", "18", "21")  # IL, IA, WI, MO, IN, KY
state_names <- c("Illinois", "Iowa", "Wisconsin", "Missouri", "Indiana", "Kentucky")

counties_sf <- counties(state = state_fips, year = 2020, cb = TRUE, progress_bar = FALSE)
counties_sf <- st_as_sf(counties_sf)

# Add 5-digit fips and state name
counties_sf$county_fips <- paste0(counties_sf$STATEFP, counties_sf$COUNTYFP)
counties_sf$state_name <- factor(counties_sf$STATEFP,
                                  levels = state_fips,
                                  labels = state_names)

cat(sprintf("Loaded %d counties across 6 states\n", nrow(counties_sf)))

# Pull state polygons for thick border lines
states_sf <- states(year = 2020, cb = TRUE, progress_bar = FALSE)
states_sf <- st_as_sf(states_sf)
states_sf <- states_sf[states_sf$STATEFP %in% state_fips, ]

# ---- 2. Load IL dispensary data -------------------------------------------

# Try common paths for your IL cross-section
xs_paths <- c(
  "data/processed/il_cross_section.rds",
  "data/processed/cross_section.rds",
  "data/processed/ch4_cross_section.rds",
  "data/processed/il_counties.rds"
)
xs_path <- xs_paths[file.exists(xs_paths)][1]

il_data <- NULL
if (!is.na(xs_path)) {
  il_data <- as.data.table(readRDS(xs_path))
  cat(sprintf("Loaded IL cross-section from %s (%d rows)\n", xs_path, nrow(il_data)))
  
  # Ensure county_fips is 5-digit character
  if ("county_fips" %in% names(il_data)) {
    il_data[, county_fips := sprintf("%05d", as.integer(county_fips))]
  } else if ("fips" %in% names(il_data)) {
    il_data[, county_fips := sprintf("%05d", as.integer(fips))]
  } else if ("FIPS" %in% names(il_data)) {
    il_data[, county_fips := sprintf("%05d", as.integer(FIPS))]
  }
} else {
  cat("WARNING: IL cross-section not found. Maps will show geometry only.\n")
  cat("Searched:", paste(xs_paths, collapse = ", "), "\n")
}

# Merge IL dispensary data onto IL counties only
il_counties <- counties_sf[counties_sf$state_name == "Illinois", ]
if (!is.null(il_data) && "EntryPer100k" %in% names(il_data)) {
  il_counties <- merge(il_counties, il_data[, .(county_fips, EntryPer100k, Border15)],
                       by = "county_fips", all.x = TRUE)
} else {
  il_counties$EntryPer100k <- NA_real_
  il_counties$Border15 <- 0
}

cat(sprintf("IL counties: %d (with dispensary data: %d)\n",
            nrow(il_counties), sum(!is.na(il_counties$EntryPer100k))))

# ---- 3. Identify Iowa treated counties (from panel) -----------------------

panel_path <- "data/processed/panel_county_year.rds"
ia_treated_fips <- c()

if (file.exists(panel_path)) {
  panel <- as.data.table(readRDS(panel_path))
  state_col <- if ("state_name" %in% names(panel)) "state_name" else "state"
  ia_treated <- unique(panel[get(state_col) %in% c("iowa", "Iowa", "IA", "ia") & border == 1,
                              county_fips])
  ia_treated_fips <- sprintf("%05d", as.integer(ia_treated))
  cat(sprintf("Iowa treated counties: %d\n", length(ia_treated_fips)))
} else {
  # Fallback: hard-code the 8 known treated Iowa Mississippi-River counties
  ia_treated_fips <- c("19045", "19057", "19061", "19097",
                        "19111", "19115", "19139", "19163")
  cat("Using hard-coded Iowa treated set (8 counties)\n")
}

# Mark Iowa counties as treated vs control
ia_counties <- counties_sf[counties_sf$state_name == "Iowa", ]
ia_counties$treated <- ia_counties$county_fips %in% ia_treated_fips

# ============================================================================
# MAP A: All 5 prohibition neighbors + Illinois
# ============================================================================

cat("\nBuilding Map A (all 5 prohibition borders)...\n")

# Set color palettes
INK <- "#1A1A1A"
MUTE <- "#777777"
ACCENT <- "#B91C1C"

map_A <- ggplot() +
  # Non-IL states in soft grey background
  geom_sf(data = counties_sf[counties_sf$state_name != "Illinois", ],
          fill = "#F4F4F4", color = "#CCCCCC", linewidth = 0.15) +
  
  # IL counties shaded by dispensary density
  geom_sf(data = il_counties,
          aes(fill = EntryPer100k), color = "white", linewidth = 0.2) +
  
  # State borders (thick)
  geom_sf(data = states_sf, fill = NA, color = INK, linewidth = 0.7) +
  
  # IL state outline emphasized
  geom_sf(data = states_sf[states_sf$STATEFP == "17", ],
          fill = NA, color = INK, linewidth = 1.1) +
  
  scale_fill_viridis(
    name = "Dispensaries\nper 100k",
    option = "magma",
    direction = -1,
    na.value = "#EEEEEE",
    breaks = c(0, 5, 10, 15, 20, 25),
    limits = c(0, NA)
  ) +
  
  # State name labels
  annotate("text", x = -89.5, y = 40.0, label = "ILLINOIS",
           fontface = "bold", size = 4.5, color = INK) +
  annotate("text", x = -93.5, y = 42.0, label = "IOWA",
           fontface = "italic", size = 3.5, color = MUTE) +
  annotate("text", x = -89.5, y = 44.5, label = "WISCONSIN",
           fontface = "italic", size = 3.5, color = MUTE) +
  annotate("text", x = -93.0, y = 38.5, label = "MISSOURI",
           fontface = "italic", size = 3.5, color = MUTE) +
  annotate("text", x = -86.0, y = 40.0, label = "INDIANA",
           fontface = "italic", size = 3.5, color = MUTE) +
  annotate("text", x = -85.0, y = 37.5, label = "KENTUCKY",
           fontface = "italic", size = 3.5, color = MUTE) +
  
  coord_sf(xlim = c(-96, -82), ylim = c(36, 47), expand = FALSE) +
  
  labs(
    title    = "Illinois dispensary density and prohibition-state neighbors",
    subtitle = "Adult-use dispensaries per 100,000 county residents (IDFPR registry, 2024)",
    caption  = "Counties: Census TIGER/Line 2020. Neighbors all prohibit recreational cannabis through 2024."
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = INK),
    plot.subtitle    = element_text(color = MUTE, size = 11),
    plot.caption     = element_text(color = MUTE, size = 8, hjust = 0),
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10, color = INK),
    legend.text      = element_text(size = 9, color = MUTE),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(1.0, "cm")
  )

ggsave("out/map_A_all_borders.pdf", map_A, width = 11, height = 8.5, device = cairo_pdf)
ggsave("out/map_A_all_borders.png", map_A, width = 11, height = 8.5, dpi = 300)

cat("Saved: out/map_A_all_borders.pdf + .png\n")

# ============================================================================
# MAP B: Iowa + Illinois (Iowa-focused)
# ============================================================================

cat("\nBuilding Map B (Iowa focus)...\n")

map_B <- ggplot() +
  # IL counties shaded by dispensary density
  geom_sf(data = il_counties,
          aes(fill = EntryPer100k), color = "white", linewidth = 0.2) +
  
  # IA control counties in light grey
  geom_sf(data = ia_counties[!ia_counties$treated, ],
          fill = "#F4F4F4", color = "#BBBBBB", linewidth = 0.15) +
  
  # IA treated counties in accent red
  geom_sf(data = ia_counties[ia_counties$treated, ],
          fill = ACCENT, color = INK, linewidth = 0.4, alpha = 0.75) +
  
  # State borders (thick black)
  geom_sf(data = states_sf[states_sf$STATEFP %in% c("17", "19"), ],
          fill = NA, color = INK, linewidth = 1.0) +
  
  scale_fill_viridis(
    name = "IL\ndispensaries\nper 100k",
    option = "magma",
    direction = -1,
    na.value = "#EEEEEE",
    breaks = c(0, 5, 10, 15, 20),
    limits = c(0, NA)
  ) +
  
  # Labels
  annotate("text", x = -89.5, y = 40.0, label = "ILLINOIS",
           fontface = "bold", size = 4.5, color = INK) +
  annotate("text", x = -93.5, y = 42.0, label = "IOWA",
           fontface = "bold", size = 4.5, color = INK) +
  annotate("text", x = -91.3, y = 41.9, label = "Mississippi River boundary",
           fontface = "italic", size = 3.0, color = ACCENT, angle = 75) +
  
  # Iowa treated counties caption
  annotate("text", x = -95.5, y = 39.5,
           label = "Red = 8 IA treated counties\n(Clinton, Des Moines, Dubuque,\nJackson, Lee, Louisa,\nMuscatine, Scott)",
           hjust = 0, size = 3.0, color = ACCENT, fontface = "italic") +
  
  coord_sf(xlim = c(-97, -86), ylim = c(36.8, 44), expand = FALSE) +
  
  labs(
    title    = "Illinois dispensary density and Iowa treated counties",
    subtitle = "Eight Iowa Mississippi-River counties form the Chapter 5 treatment group",
    caption  = "Counties: Census TIGER/Line 2020.  IL dispensaries: IDFPR registry, 2024."
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, color = INK),
    plot.subtitle    = element_text(color = MUTE, size = 11),
    plot.caption     = element_text(color = MUTE, size = 8, hjust = 0),
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 10, color = INK),
    legend.text      = element_text(size = 9, color = MUTE),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(1.0, "cm")
  )

ggsave("out/map_B_iowa_focus.pdf", map_B, width = 11, height = 8, device = cairo_pdf)
ggsave("out/map_B_iowa_focus.png", map_B, width = 11, height = 8, dpi = 300)

cat("Saved: out/map_B_iowa_focus.pdf + .png\n")

# ============================================================================
# Sanity checks
# ============================================================================

cat("\n=== Sanity checks ===\n")
cat(sprintf("IL counties total:                %d\n", nrow(il_counties)))
cat(sprintf("IL counties with positive entry:  %d\n",
            sum(il_counties$EntryPer100k > 0, na.rm = TRUE)))
cat(sprintf("IL Border15 counties (if data):   %d\n",
            sum(il_counties$Border15 == 1, na.rm = TRUE)))
cat(sprintf("IA counties total:                %d\n", nrow(ia_counties)))
cat(sprintf("IA treated (Mississippi border):  %d\n", sum(ia_counties$treated)))

cat("\nDone. Use map_A for full-thesis context; map_B for Ch5 focus.\n")
cat("Both rendered at PDF + PNG (300 DPI).\n")
