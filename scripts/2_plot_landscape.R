library(gridExtra)
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions") - only needed to save data
library(raster)
library(rasterVis)
library(sp)
library(tidyverse)

# Import sampling points
sampling_points <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))
names(sampling_points) <- c("x", "y", "Group", "Site_name") # rename columns

# Import NLCD layer
nlcd <- readRDS(paste0(getwd(), "/data/output/nlcd.rds"))
nlcd_reclassified <- readRDS(paste0(getwd(), "/data/output/nlcd_reclassified.rds"))

# Import habitat surface
habitat_surface <- raster::raster(paste0(getwd(), "/data/GIS/habitat_surface.tif"))

# plot original NLCD
plot_nlcd <- rasterVis::levelplot(nlcd, 
                                  main = "NLCD") + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="green", size = 2))

# plot reclassified NLCD
plot_nlcd_reclassified <- rasterVis::levelplot(nlcd_reclassified, 
                                               main = "NLCD reclassified") + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="green", size = 2))

# plot original vs reclassified NLCD
plot_nlcd_overall <- gridExtra::grid.arrange(plot_nlcd, 
                                             plot_nlcd_reclassified, ncol = 2)


# plot habitat surface
plot_surface <- rasterVis::levelplot(habitat_surface, main = "Surface metrics", 
                                     margin = FALSE) + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="green", size = 2))

# plot habitat surface vs. reclassified NLCD``
plot_surface_patch <- gridExtra::grid.arrange(plot_nlcd_reclassified, 
                                              plot_surface, ncol = 2)
