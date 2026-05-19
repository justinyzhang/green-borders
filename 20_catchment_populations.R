# =============================================================================
# 20_catchment_populations.R
# =============================================================================
# Foundational GIS preprocessing for the radius design (Lychagin 5/14 directive).
# For each focal county, computes total population within radius r in miles,
# partitioned by source-state group.
#
# Two focal sets:
#   (a) IL counties (102) for the Chapter 4 gravity entry regression.
#   (b) IA counties (99) for the Chapter 5 continuous-dose DiD on NIBRS arrests.
#
# Three radii: 15, 50, 100 miles, corresponding to immediate-border interaction,
# one-county catchment, and two-county gravity catchment respectively (the
# average IL county has linear extent ~13 miles, so 15/50/100 maps to
# within-county / adjacent ring / two-ring catchments).
#
# Source: Census Centers of Population 2020 (county pop-weighted centroids).
# Population: ACS 5-year 2018-2022 block-group totals.
# CRS: EPSG 5070 Conterminous US Albers (meters; preserves area, true distance).
#
# Outputs:
#   out/catchment_pop_il_focal.csv  (306 rows: 102 counties x 3 radii)
#   out/catchment_pop_ia_focal.csv  (297 rows: 99 counties x 3 radii)
#
# Runtime: ~15-25 minutes depending on machine; sf operations on ~40k block
# groups are the bottleneck.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(tigris)
  library(tidycensus)
})

options(tigris_use_cache = TRUE, tigris_class = "sf")
sf::sf_use_s2(FALSE)  # Use planar geometry on projected CRS for buffer distance accuracy

# ---- 0. Configuration -------------------------------------------------------

TARGET_CRS <- 5070               # Conterminous US Albers Equal Area Conic (meters)
MILES_TO_METERS <- 1609.344
ACS_YEAR <- 2022                 # ACS 5-year 2018-2022
RADII_MILES <- c(15, 50, 100)
OUT_DIR <- "out"
dir.create(OUT_DIR, showWarnings = FALSE)

# State sets relevant to the radius design.
# Source states: anything within 100 miles of any IL or IA county centroid.
STATES_RELEVANT <- c("IL","IA","MO","IN","KY","WI","MI","TN","MN","NE","KS","OK")

# State-level prohibition status as of January 2020 (the IL treatment shock).
# Note: MI legalized 2018 but does not share a land border with IL.
STATE_FIPS <- tribble(
  ~abbr, ~fips,  ~status_2020,
  "IL",  "17",   "legal",
  "IA",  "19",   "prohibition",
  "MO",  "29",   "prohibition",   # legalized Feb 2023; held out of primary window
  "IN",  "18",   "prohibition",
  "KY",  "21",   "prohibition",
  "WI",  "55",   "prohibition",
  "MI",  "26",   "legal",
  "TN",  "47",   "prohibition",
  "MN",  "27",   "prohibition",   # legalized Aug 2023
  "NE",  "31",   "prohibition",
  "KS",  "20",   "prohibition",
  "OK",  "40",   "prohibition"
)

# ---- 1. Pull block-group populations from ACS ------------------------------

message("[1] Loading ACS 2018-2022 5-year block-group populations for ",
        length(STATES_RELEVANT), " states. This may take a few minutes.")

# census_api_key("YOUR_KEY_HERE", install = TRUE) # Run this once if not set

bg_pop_list <- vector("list", length(STATES_RELEVANT))
names(bg_pop_list) <- STATES_RELEVANT
for (st in STATES_RELEVANT) {
  message("    Pulling ", st, " ...")
  bg_pop_list[[st]] <- get_acs(
    geography = "block group",
    variables = "B01003_001",        # total population estimate
    state = st,
    year = ACS_YEAR,
    survey = "acs5",
    geometry = TRUE,
    output = "wide",
    cache_table = TRUE
  )
}

bg_pop <- bind_rows(bg_pop_list) %>%
  rename(pop = B01003_001E) %>%
  select(GEOID, pop, geometry) %>%
  mutate(state_fips = substr(GEOID, 1, 2)) %>%
  st_transform(TARGET_CRS)

message("    Loaded ", nrow(bg_pop), " block groups.")

# Block-group centroid points: faster intersection check than polygons.
# Population-weighted centroid is unavailable at BG level, so we use geographic
# centroid; the error is small relative to the radius scales (15-100 miles).
bg_centroids <- bg_pop %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  select(GEOID, pop, state_fips)

# ---- 2. Load county population-weighted centroids --------------------------

message("[2] Loading Census 2020 population-weighted county centroids.")

cop_url <- "https://www2.census.gov/geo/docs/reference/cenpop2020/county/CenPop2020_Mean_CO.txt"
cop <- read_csv(
  cop_url,
  col_types = cols(
    STATEFP = "c", COUNTYFP = "c", COUNAME = "c", STNAME = "c",
    POPULATION = "d", LATITUDE = "d", LONGITUDE = "d"
  )
) %>%
  mutate(GEOID = paste0(STATEFP, COUNTYFP))

cop_sf <- cop %>%
  st_as_sf(coords = c("LONGITUDE","LATITUDE"), crs = 4326) %>%
  st_transform(TARGET_CRS)

# ---- 3. Catchment function -------------------------------------------------

#' Compute catchment population for each focal county.
#'
#' @param focal_centroids sf POINT object with GEOID column
#' @param bg sf POINT object (block-group centroids) with GEOID, pop, state_fips
#' @param radii numeric vector of radii in miles
#' @param state_groups named list of state_fips vectors, e.g.
#'   list(IL = "17", PROH = c("19","29","18","21","55"))
#' @return tibble with columns: GEOID, radius_mi, pop_<group> for each state group
compute_catchment <- function(focal_centroids, bg, radii, state_groups) {

  result <- focal_centroids %>%
    st_drop_geometry() %>%
    select(GEOID) %>%
    expand_grid(radius_mi = radii)

  for (sg_name in names(state_groups)) {
    sg_fips <- state_groups[[sg_name]]
    bg_sg <- bg %>% filter(state_fips %in% sg_fips)
    message("    Group '", sg_name, "': ", nrow(bg_sg), " block groups; ",
            "states ", paste(sg_fips, collapse = ","))

    # For each radius, find which bg centroids are within distance of each focal
    radius_results <- map_dfr(radii, function(r) {
      dist_m <- r * MILES_TO_METERS
      # st_is_within_distance returns a list with one element per focal county
      within_idx <- st_is_within_distance(focal_centroids, bg_sg, dist = dist_m)

      tibble(
        GEOID = focal_centroids$GEOID,
        radius_mi = r,
        pop_sum = vapply(within_idx, function(idx) {
          if (length(idx) == 0) 0 else sum(bg_sg$pop[idx], na.rm = TRUE)
        }, numeric(1))
      )
    })

    col_name <- paste0("pop_", sg_name)
    result <- result %>%
      left_join(radius_results %>% rename(!!col_name := pop_sum),
                by = c("GEOID","radius_mi"))
  }

  result
}

# ---- 4. Define focal sets and source state groups --------------------------

# 4a. IL focal: all 102 IL counties (entry regression sample).
il_focal <- cop_sf %>% filter(STATEFP == "17")
stopifnot(nrow(il_focal) == 102)

# 4b. IA focal: all 99 IA counties (crime DiD sample).
ia_focal <- cop_sf %>% filter(STATEFP == "19")
stopifnot(nrow(ia_focal) == 99)

# 4c. Cross-state placebo focal (for Chapter 4 column 3):
#     Counties in non-IL prohibition states that border OTHER prohibition states.
#     Example: IA counties bordering MO (pre-2023), IN counties bordering KY, etc.
#     We pull all counties in IA, MO, IN, KY, WI and let downstream code restrict.
placebo_focal <- cop_sf %>% filter(STATEFP %in% c("19","29","18","21","55"))

# Source-state groups for the IL focal (entry regression):
# - IL-side population (within-state market)
# - Prohibition-side population (cross-border market; IL's 5 prohibition land
#   neighbors: IA, MO, IN, KY, WI)
src_groups_il <- list(
  IL   = STATE_FIPS$fips[STATE_FIPS$abbr == "IL"],
  PROH = STATE_FIPS$fips[STATE_FIPS$abbr %in% c("IA","MO","IN","KY","WI")]
)

# Source-state groups for the IA focal (crime DiD):
# - IL-side population (legal-state dose; the treatment intensity)
# - IA-side population (within-state baseline catchment)
src_groups_ia <- list(
  IL = STATE_FIPS$fips[STATE_FIPS$abbr == "IL"],
  IA = STATE_FIPS$fips[STATE_FIPS$abbr == "IA"]
)

# Source-state groups for the placebo (cross-state):
# For each placebo focal county, we want population within radius of "other
# prohibition states" (not its own state). The grouping logic is handled
# downstream by joining each focal county to its neighboring states' pop.
src_groups_placebo <- list(
  PROH_NEIGHBORS = STATE_FIPS$fips[STATE_FIPS$status_2020 == "prohibition"]
)

# ---- 5. Compute catchments -------------------------------------------------

message("[3] Computing catchments for IL focal (102 counties, 3 radii) ...")
cat_il <- compute_catchment(il_focal, bg_centroids, RADII_MILES, src_groups_il) %>%
  left_join(cop %>% select(GEOID, COUNAME, STNAME, county_pop_2020 = POPULATION),
            by = "GEOID") %>%
  select(GEOID, COUNAME, STNAME, county_pop_2020, radius_mi, pop_IL, pop_PROH)

message("[4] Computing catchments for IA focal (99 counties, 3 radii) ...")
cat_ia <- compute_catchment(ia_focal, bg_centroids, RADII_MILES, src_groups_ia) %>%
  left_join(cop %>% select(GEOID, COUNAME, STNAME, county_pop_2020 = POPULATION),
            by = "GEOID") %>%
  select(GEOID, COUNAME, STNAME, county_pop_2020, radius_mi, pop_IL, pop_IA)

message("[5] Computing catchments for placebo focal (cross-state) ...")
cat_placebo <- compute_catchment(placebo_focal, bg_centroids, RADII_MILES,
                                 src_groups_placebo) %>%
  left_join(cop %>% select(GEOID, COUNAME, STNAME, county_pop_2020 = POPULATION),
            by = "GEOID") %>%
  select(GEOID, COUNAME, STNAME, county_pop_2020, radius_mi, pop_PROH_NEIGHBORS)

# ---- 6. Write outputs ------------------------------------------------------

write_csv(cat_il,      file.path(OUT_DIR, "catchment_pop_il_focal.csv"))
write_csv(cat_ia,      file.path(OUT_DIR, "catchment_pop_ia_focal.csv"))
write_csv(cat_placebo, file.path(OUT_DIR, "catchment_pop_placebo_focal.csv"))

message("\n[DONE] Outputs written:")
message("  ", file.path(OUT_DIR, "catchment_pop_il_focal.csv"),
        " (", nrow(cat_il), " rows)")
message("  ", file.path(OUT_DIR, "catchment_pop_ia_focal.csv"),
        " (", nrow(cat_ia), " rows)")
message("  ", file.path(OUT_DIR, "catchment_pop_placebo_focal.csv"),
        " (", nrow(cat_placebo), " rows)")

# Quick sanity checks
cat("\n--- Sanity check: IL focal at r = 15 mi ---\n")
print(cat_il %>% filter(radius_mi == 15) %>%
        arrange(desc(pop_PROH)) %>% head(10))

cat("\n--- Sanity check: IA focal at r = 15 mi ---\n")
print(cat_ia %>% filter(radius_mi == 15) %>%
        arrange(desc(pop_IL)) %>% head(10))

cat("\n--- Sanity check: IA bridge counties (Scott, Dubuque, etc.) ---\n")
bridge_geoids <- c("19045","19057","19061","19097","19111","19139","19163","19115")
print(cat_ia %>% filter(GEOID %in% bridge_geoids, radius_mi == 50) %>%
        arrange(desc(pop_IL)))
