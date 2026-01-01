## 0. Setup -------------------------------------------------------------------

# install.packages("here")   # Install once per machine
library(here)

library(terra)
library(sf)
library(dplyr)
library(tibble) 
library(readr)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------#
#
# This repository includes a *cropped GeoTIFF derived from the Copernicus/ESA
# C3S global land-cover product (NetCDF-4; 300 m), using only the "lccs_class"
# layer for year 2022 (product version 2.1.1).
# - The full global NetCDF is very large and is not supported in this repo.
# - For this analysis we only need land cover around the site midpoints, so we
#   provide a pre-cropped raster that covers the union of 5 km site buffers (plus a
#   500 m margin), generated once offline and committed to GitHub via Git LFS.
#
# To access the full global product if required, download the
# official dataset from the Copernicus Climate Data Store (CDS) and use the
# NetCDF "lccs_class" variable directly (see README documentation).
# -----------------------------------------------------------------------------#

lc_file <- here::here(
  "Dataframes", "Land use mapping",
  "lccs_class_2022_v2.1.1_crop_sites_5kmplus500m.tif"
)

# Load raster (already 'lccs_class')
lc_raster <- terra::rast(lc_file)
print(lc_raster)

## 1. Define 'human-dominated' and 'water' integer codes -----------------------

# Following the ESA CCI LCCS legend:
#   10–40   = Cropland subclasses
#   190     = Urban / built‐up
#   210     = Permanent water bodies
#   220     = Coastal water (ocean)
human_classes <- c(10:40, 190)
water_classes <- c(210, 220)

# Mask out water pixels so that only land codes remain
#   Both stages will take some time
water_mask <- lc_raster %in% water_classes

lc_terrestrial <- mask(
  lc_raster,
  water_mask,
  maskvalues  = TRUE,
  updatevalue = NA
)

## 2. Process midpoint coordinates ---------------------------------------------

# Load midpoints
coords_df <- tribble(
  ~region      ,  ~coast_lat   ,  ~coast_lon,
  "Argentina"  ,  -54.83997   ,  -68.323227,
  "Chile"      ,  -41.49617   ,  -72.985483,
  "Denmark"    ,   55.08500   ,    8.568000,
  "Ecuador"    ,   -0.74337   ,  -90.303910,
  "Greenland"  ,   64.17800   ,  -51.701000,
  "New Zealand",  -41.23803   ,  174.897015,
  "Norway"     ,   69.64420   ,   18.951990,
  "Portugal1"  ,   32.64543   ,  -16.911290,
  "Portugal2"  ,   41.56397   ,   -8.798142,
  "UK"         ,   50.78500   ,    0.021000,
  "Vietnam1"   ,   20.72611   ,  107.045390
)

sites_sf <- coords_df %>%
  st_as_sf(coords = c("coast_lon", "coast_lat"), crs = 4326)

# Build 5 km buffers around each site (requires a metre‐based CRS)
if (is.lonlat(lc_terrestrial)) {
  # If raster is in geographic (degrees); reproject sites to EPSG:3857 for buffering
  sites_proj    <- st_transform(sites_sf, crs = 3857)   # Web Mercator (units = meters)
  buffers_5km_m <- st_buffer(sites_proj, dist = 5000)   # 5000 m
  # Reproject buffers back to raster CRS (EPSG:4326)
  buffers_5km   <- st_transform(buffers_5km_m, crs = crs(lc_terrestrial))
} else {
  # If lc_terrestrial were already in a projected CRS (unlikely), buffer directly
  sites_proj  <- st_transform(sites_sf, crs = crs(lc_terrestrial))
  buffers_5km <- st_buffer(sites_proj, dist = 5000)
}

# Crop lc_terrestrial to only the area around all buffers (±500 m margin)
# Compute all buffers, then add a 500 m margin
union_bbox <- st_union(buffers_5km)
bbox_small <- st_buffer(union_bbox, dist = 500)

# Convert bbox_small to SpatExtent by first converting to SpatVector
bbox_small_proj <- st_transform(bbox_small, crs = crs(lc_terrestrial))
bbox_vect       <- vect(bbox_small_proj)               # Convert to SpatVector
terra_extent    <- ext(bbox_vect)                      # SpatExtent for cropping

# Crop lc_terrestrial by the small extent (will take a bit of time)
lc_crop_small <- crop(lc_terrestrial, terra_extent)

# Inspect to ensure lc_crop_small is smaller than lc_terrestrial
print(lc_terrestrial)
print(lc_crop_small)

###############################################################################

## 3. Extract % human‐dominated land within each 5 km buffer --------------------

# Convert the 5 km buffers (sf) to a SpatVector in the same CRS as lc_crop_small
buffers_vect <- vect(st_transform(buffers_5km, crs = crs(lc_crop_small)))

# Use terra::extract() to count only land pixels (non‐NA)
pct_human_df <- terra::extract(
  lc_crop_small,
  buffers_vect,
  fun = function(vals) {
    total_land_pixels <- sum(!is.na(vals))
    if (total_land_pixels == 0) {
      return(NA_real_)  # No land pixels under this buffer
    }
    human_pixels <- sum(vals %in% human_classes, na.rm = TRUE)
    (human_pixels / total_land_pixels) * 100
  }
)

# Rename new column and merge with the original site IDs
pct_col_name <- names(pct_human_df)[2]   # usually “lyr.1”
pct_human_df <- pct_human_df %>%
  rename(pct_human_landcover_5km = all_of(pct_col_name))

site_ids_df <- sites_sf %>%
  st_drop_geometry() %>%
  mutate(ID = seq_len(n()))

result_table <- site_ids_df %>%
  left_join(pct_human_df, by = "ID") %>%
  dplyr::select(region, pct_human_landcover_5km)

print(result_table)

## 4. Mapping -------------------------------------------------------------------

# Define six land‐cover categories and colors
category_colors <- c(
  "Cropland"      = "#D73027",
  "Urban"         = "#FC8D59",
  "Forest"        = "#1A9850",
  "Shrub/Grass"   = "#91CF60",
  "Sparse/Barren" = "#E0E0E0",
  "Other"         = "#999999"
)
all_categories <- names(category_colors)

# Process site labels so 'Portugal1' becomes 'Portugal 1' etc.
sites_sf$label <- sub("([A-Za-z]+)([0-9]+)$", "\\1 \\2", sites_sf$region)

plot_list <- vector("list", length = nrow(sites_sf))

for (i in seq_len(nrow(sites_sf))) {
  site_label <- sites_sf$label[i]
  site_pt    <- sites_sf[i, ]
  
  # Build the 5 km buffer (EPSG:3857 → EPSG:4326)
  site_3857    <- st_transform(site_pt, crs = 3857)
  buffer_3857  <- st_buffer(site_3857, dist = 5000)
  buffer_wgs84 <- st_transform(buffer_3857, crs = 4326)
  
  # Build a 500 m margin in degrees for cropping
  buffer_margin <- st_buffer(buffer_wgs84, dist = 0.005)
  
  # Crop lc_crop_small to this margin
  margin_vect <- vect(st_transform(buffer_margin, crs = crs(lc_crop_small)))
  terra_ext2  <- ext(margin_vect)
  lc_crop2    <- crop(lc_crop_small, terra_ext2)
  
  # Mask lc_crop2 by the 5 km circle
  lc_crop_masked <- mask(
    lc_crop2,
    vect(st_transform(buffer_wgs84, crs = crs(lc_crop2)))
  )
  
  # Convert masked raster to a data.frame, dropping NA
  df_lc_inside <- as.data.frame(
    lc_crop_masked,
    xy    = TRUE,
    na.rm = TRUE
  )
  names(df_lc_inside)[3] <- "class"
  
  # Map integer code to one of six categories
  df_lc_inside$category <- factor(
    case_when(
      df_lc_inside$class %in% 10:40                  ~ "Cropland",
      df_lc_inside$class == 190                      ~ "Urban",
      df_lc_inside$class %in% c(50, 60, 70, 100, 110) ~ "Forest",
      df_lc_inside$class %in% c(121, 130, 140, 150)  ~ "Shrub/Grass",
      df_lc_inside$class %in% c(160, 170)            ~ "Sparse/Barren",
      TRUE                                           ~ "Other"
    ),
    levels = all_categories
  )
  
  # Build subplot
  p <- ggplot() +
    # Blue (ocean) fill for entire 5 km circle
    geom_sf(
      data        = buffer_wgs84,
      inherit.aes = FALSE,
      fill        = "lightblue",
      color       = NA
    ) +
    # Coloured land pixels above the blue
    geom_raster(
      data = df_lc_inside,
      aes(x = x, y = y, fill = category)
    ) +
    # Black outline of the 5 km circle
    geom_sf(
      data        = buffer_wgs84,
      inherit.aes = FALSE,
      fill        = NA,
      color       = "black",
      size        = 1
    ) +
    # White dot = site midpoint
    geom_sf(
      data        = site_pt,
      inherit.aes = FALSE,
      shape       = 21,
      fill        = "white",
      color       = "black",
      size        = 3,
      stroke      = 0.8
    ) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    coord_sf(
      xlim   = st_bbox(buffer_margin)[c("xmin", "xmax")],
      ylim   = st_bbox(buffer_margin)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    labs(title = site_label, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title      = element_text(size = 10, face = "bold"),
      legend.position = "none",
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank()
    )
  
  plot_list[[i]] <- p
}

legend_df <- data.frame(
  category = factor(all_categories, levels = all_categories),
  x = 1,
  y = seq_along(all_categories)
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, fill = category)) +
  geom_tile(width = 0, height = 0, show.legend = TRUE) +
  scale_fill_manual(
    name   = "Land cover category",
    values = category_colors,
    drop   = FALSE,
    breaks = all_categories,
    labels = all_categories
  ) +
  guides(fill = guide_legend(
    title.position = "top",
    title.hjust    = 0.5,
    ncol           = 1,
    byrow          = TRUE
  )) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 8),
    legend.key.height = unit(0.7, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    panel.background  = element_blank(),
    plot.background   = element_blank()
  )

# Combine the subplots and legend
panel <- (
  wrap_plots(
    plot_list,
    ncol          = 4,
    align         = "hv",
    widths        = rep(1, 4),
    panel_spacing = unit(0.2, "cm")
  )
  |
    legend_plot
) +
  plot_layout(widths = c(4, 0.2)) +
  plot_annotation(title = "5 km land Cover (2022) for coastal midpoints") +
  theme(
    plot.title    = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

# Inset north arrow
north_plot <- ggplot() +
  annotate(
    "segment",
    x = 0.5, xend = 0.5,
    y = 0.2, yend = 0.8,
    arrow = arrow(length = unit(0.2, "cm"))
  ) +
  annotate("text", x = 0.5, y = 0.9, label = "N", size = 6, fontface = "bold") +
  theme_void()

panel_with_north <- panel +
  inset_element(
    north_plot,
    left   = 0.78, bottom = 0.80,
    right  = 0.97, top    = 0.97,
    on_top = TRUE
  )

print(panel_with_north)