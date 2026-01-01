## 0. Setup -------------------------------------------------------------------

# install.packages("here")   # Install once per machine should code and data be read from GitHub
library(here)

library(sf)
library(dplyr)
library(raster)
library(tibble)
library(purrr)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(viridis)  
library(ggspatial)
library(patchwork)  
library(prettymapr)
library(elevatr)
library(ggpattern)

# Polluted/non-polluted coordinates
coords_df <- tribble(
  ~region,      ~site,          ~lat,         ~lon,
  "Chile",      "Polluted",    -41.496167,   -72.985483,
  "Chile",      "Non-Polluted",-41.594317,   -72.705000,
  "Argentina",  "Polluted",    -54.839968,   -68.323227,
  "Argentina",  "Non-Polluted",-54.818830,   -68.181426,
  "Ecuador",    "Polluted",     -0.743370,   -90.303910,
  "Ecuador",    "Non-Polluted",-0.765930,   -90.340570,
  "New Zealand","Polluted",    -41.238031,   174.897015,
  "New Zealand","Non-Polluted",-41.320606,   174.870620,
  "UK",         "Polluted",     50.785000,     0.021000,
  "UK",         "Non-Polluted", 50.789000,    -0.002000,
  "Denmark",    "Polluted",     55.085000,     8.568000,
  "Denmark",    "Non-Polluted", 55.146000,     8.606000,
  "Norway",     "Polluted",     69.644200,    18.951990,
  "Norway",     "Non-Polluted", 69.634570,    18.902590,
  "Portugal1",  "Polluted",     32.645434,   -16.911290,
  "Portugal1",  "Non-Polluted", 32.741415,   -16.714039,
  "Portugal2",  "Polluted",     41.563968,    -8.798142,
  "Portugal2",  "Non-Polluted", 41.309617,    -8.742615,
  "Greenland",  "Polluted",     64.178000,   -51.701000,
  "Greenland",  "Non-Polluted", 64.167000,   -51.746000,
  "Vietnam1",   "Polluted",     20.726110,   107.045390,
  "Vietnam1",   "Non-Polluted", 20.740210,   106.997420
)

# convert to sf
sites_sf <- coords_df %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326)

## 1. Compute true geodesic midpoint for each region -------------------------------------------------------------------

mid_df <- coords_df %>%
  group_by(region) %>%
  summarize(
    lon1 = lon[site=="Polluted"],
    lat1 = lat[site=="Polluted"],
    lon2 = lon[site=="Non-Polluted"],
    lat2 = lat[site=="Non-Polluted"],
    .groups="drop"
  ) %>%
  rowwise() %>%
  mutate(
    mid = midPoint(c(lon1, lat1), c(lon2, lat2)),
    gc_lon = mid[1],
    gc_lat = mid[2]
  ) %>%
  ungroup() %>%
  dplyr::select(region, gc_lon, gc_lat)

## 2. Load global R coastline -------------------------------------------------------------------

coast_global <- ne_coastline(scale="medium", returnclass="sf")

## 3. Align midpoint to its local mainland coast -------------------------------------------------------------------

snap_to_coast <- function(region_name, lon, lat) {
  # Original two sites as sf
  pts <- coords_df %>%
    filter(region == region_name) %>%
    st_as_sf(coords = c("lon","lat"), crs = 4326)
  
  # Small bounding box and buffer around them
  local_bb    <- st_bbox(pts) %>% st_as_sfc() %>% st_buffer(0.5)
  local_coast <- st_crop(coast_global, local_bb)
  if (nrow(local_coast)==0) local_coast <- coast_global
  
  # Midpoint as sf point
  mid_pt_sfc  <- st_sfc(st_point(c(lon, lat)), crs = 4326)
  
  # Select the nearest coastline feature
  nearest_idx <- st_nearest_feature(mid_pt_sfc, local_coast)
  one_coast   <- local_coast[nearest_idx, ]
  
  # Get the line from midpoint to the feature, then its endpoint
  nbr_line    <- st_nearest_points(mid_pt_sfc, one_coast)
  coast_pt    <- st_cast(nbr_line, "POINT")[2]
  
  # Extract extract coordinates
  xy <- st_coordinates(coast_pt)
  tibble(
    region    = region_name,
    coast_lon = xy[1],
    coast_lat = xy[2]
  )
}

coastal_midpoints <- mid_df %>%
  dplyr::select(region, gc_lon, gc_lat) %>%
  purrr::pmap_dfr(~ snap_to_coast(..1, ..2, ..3))

print(coastal_midpoints)

## 4. Load climate classifications -------------------------------------------------------------------

# Load 0.01° GeoTIFF (1991–2020) (Beck et al. 2023) https://www.nature.com/articles/s41597-023-02549-6
base_dir <- here::here("Dataframes", "Land use mapping", "kg_tif", "1991_2020")
tif_file <- list.files(base_dir, pattern="koppen_geiger_0p01\\.tif$", full.names=TRUE)[1]
kg_r     <- raster(tif_file)

# Ratify and assign the 30 classifications
kg_r <- ratify(kg_r)
rat  <- levels(kg_r)[[1]]

legend_codes <- c(
  "Af",  # 1 Tropical rainforest
  "Am",  # 2 Tropical monsoon
  "Aw",  # 3 Tropical savannah
  "BWh", # 4 Hot desert
  "BWk", # 5 Cold desert
  "BSh", # 6 Hot steppe
  "BSk", # 7 Cold steppe
  "Csa", # 8 Med. dry-hot summer
  "Csb", # 9 Med. dry-warm summer
  "Csc", # 10 Med. dry-cold summer
  "Cwa", # 11 Monsoon-influenced humid subtrop., dry winter, hot summer
  "Cwb", # 12 Monsoon-influenced humid subtrop., dry winter, warm summer
  "Cwc", # 13 Monsoon-influenced humid subtrop., dry winter, cool summer
  "Cfa", # 14 Humid subtrop., no dry season, hot summer
  "Cfb", # 15 Oceanic, no dry season, warm summer
  "Cfc", # 16 Subpolar oceanic, no dry season, cool summer
  "Dsa", # 17 Cold, dry summer, hot summer
  "Dsb", # 18 Cold, dry summer, warm summer
  "Dsc", # 19 Cold, dry summer, cold summer
  "Dsd", # 20 Cold, dry summer, very cold winter
  "Dwa", # 21 Cold, dry winter, hot summer
  "Dwb", # 22 Cold, dry winter, warm summer
  "Dwc", # 23 Cold, dry winter, cold summer
  "Dwd", # 24 Cold, dry winter, very cold winter
  "Dfa", # 25 Cold, no dry season, hot summer
  "Dfb", # 26 Cold, no dry season, warm summer
  "Dfc", # 27 Cold, no dry season, cold summer
  "Dfd", # 28 Cold, no dry season, very cold winter
  "ET",  # 29 Polar tundra
  "EF"   # 30 Polar frost
)

# Check
if (length(legend_codes) != nrow(rat)) {
  stop("Your legend_codes length (", length(legend_codes),
       ") does not match your raster’s IDs (", nrow(rat), ").")
}

rat$climate    <- legend_codes[ rat$ID ]
levels(kg_r)   <- rat

## 5. Extract classification at each coastal midpoint -------------------------------------------------------------------
pts_mat <- as.matrix(coastal_midpoints[, c("coast_lon","coast_lat")])
codes <- raster::extract(kg_r, pts_mat)

# For those that fall exactly on ocean pixels, assign the modal within 10 km
na_i    <- which(is.na(codes))
if (length(na_i)) {
  buf_vals <- raster::extract(kg_r,
                      pts_mat[na_i, , drop=FALSE],
                      buffer = 10000,
                      fun    = modal,
                      na.rm  = TRUE)
  codes[na_i] <- sapply(buf_vals, `[`, 1)
}

# Apply
coastal_midpoints <- coastal_midpoints %>%
  mutate(
    koppen_code = codes,
    koppen      = factor(koppen_code,
                         levels = rat$ID,
                         labels = rat$climate)
  )

print(coastal_midpoints)

## 6. Mapping -------------------------------------------------------------------

# R Natural earth world country polygons
world <- ne_countries(scale="medium", returnclass="sf")

mid_sf <- st_as_sf(
  coastal_midpoints,
  coords = c("coast_lon","coast_lat"),
  crs    = 4326
)

worldplot <- ggplot(world) +
  geom_sf(fill   = "#f3f4ed", colour = "grey70") +
  geom_sf(
    data   = mid_sf,
    aes(color = koppen),
    shape  = 21,
    fill   = NA,
    size   = 4.5,
    stroke = 1.7
  ) +
  scale_color_viridis_d(name = "Köppen–Geiger\nclimate zone", na.value = "black") +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  coord_sf(expand = FALSE) +
  labs(
    title   = "Experimental locations\n(paired site midpoints)",
    caption = "Climate data & classes after Beck et al. (2023)",
    x       = "Longitude (°)",
    y       = "Latitude (°)"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 3)
  )) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.background = element_rect(fill = "#abc7cc", colour = NA),
    panel.grid.major = element_line(color = "grey60", linewidth = 0.1, linetype = "solid"),
    axis.text        = element_text(color = "grey20"),
    axis.ticks       = element_line(color = "grey20"),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size =  8),
    legend.key.size  = unit(0.5, "cm")
  )

plot(worldplot)