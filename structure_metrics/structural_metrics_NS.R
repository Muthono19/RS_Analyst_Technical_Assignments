# LiDAR Vegetation Structure Metrics
# Author: Stella Muthoni Gachoki
# Purpose: Derive vegetation structure metrics from .las file

# Load packages
#library(annotater) # Annotate Package Load Calls
library(lidR) # Airborne LiDAR Data Manipulation and Visualization for Forestry
library(terra) # Spatial Data Analysis
library(rgl) # 3D Visualization Using OpenGL
library(sf) # Simple Features for R
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(rlang) # Functions for Base Types and Core R and 'Tidyverse' Features

setwd("C:/Users/gstel/OneDrive/Desktop/My documents/Others/Natural State Assignment/Structural_metrics")

# Add the full LiDAR data
las <- readLAS("ALS_pointcloud.las")

# Extract Lidar data extent and plot
las_extent <- st_as_sfc(st_bbox(las), crs = 32637)
las_extent_sf <- st_sf(id = 1, geometry = las_extent)
#st_write(las_extent_sf, "las_extent.shp", delete_layer = TRUE)
st_area(las_extent_sf) 

# for computation purposes Slice the extent into 4 equal tiles
# and select one tile to use as clip object
tiles <- st_make_grid(las_extent_sf, n = c(2, 2), what = "polygons")
tiles_sf <- st_sf(tile_id = 1:length(tiles), geometry = tiles)
roi <- tiles_sf[1, ]  # or [2,], [3,], [4,] to pick other tiles

#Plot to verify it overlaps and identify the selected tile
plot(st_geometry(las_extent_sf), border = 'blue', lwd = 2, main = "LAS Extent Sliced into 4 Tiles")
plot(st_geometry(tiles_sf), add = TRUE, border = 'gray')
plot(st_geometry(roi), add = TRUE, border = 'red', lwd = 2)

# clip the original LAS the selected grid
las_clip <- clip_roi(las, roi)
las_clip

#Normalize heights using TIN ground model
las_norm <- normalize_height(las_clip, tin())
gc()

# calculate selected metrics that are more relevant to savanna grassland set up
metrics <- grid_metrics(las_norm, ~list(
  zmax = max(Z),  #Maximum height of vegetation within 2m pixel
  zq95 = quantile(Z, 0.95),  # robust canopy structure
  zsd = sd(Z),  #vertical heterogeneity (structural complexity)
  point_density = length(Z),  # horizontal vegetation density or "thickness"
  canopy_cover_tall = sum(Z > 2) / length(Z),  #tall trees (actual canopy)
  canopy_cover_all_no_ground = sum(Z > 0.2) / length(Z),  #All vegetation with no ground points
  low_veg_no_ground = sum(Z > 0.2 & Z <= 1.37) / length(Z)  #low vegetation excluding ground points
), res = 2)  # Grid resolution: 2 meters

# save the plot
png(file = "Structure_metrics.png", width = 4500, height =4500, units = "px", res = 700, type = "cairo")
plot(metrics)
dev.off()

###------------------ THE END ----------------------------------------------------------




