library(gridExtra)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt") - only needed to save data
library(raster)
library(rasterVis)
library(sp)
library(tidyverse)

# import sampling points
sampling_points <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))
names(sampling_points) <- c("x", "y", "Group", "Site_name") # rename columns

# import NLCD layer
nlcd <- readRDS(paste0(getwd(), "/data/output/nlcd.rds"))

# import color scheme
nlcd_color <- read.csv(paste0(getwd(), "/data/GIS/nlcd_color.csv"), 
                       stringsAsFactors = FALSE)

# add color scheme to RAT
levels_nlcd <- levels(nlcd)[[1]]
levels_nlcd <- dplyr::left_join(x = levels_nlcd, y = nlcd_color,
                                by = c("ID" = "value"))

levels(nlcd) <- levels_nlcd

# import reclassified NLCD layer
nlcd_reclassified <- readRDS(paste0(getwd(), "/data/output/nlcd_reclassified.rds"))

# import habitat surface
habitat_surface <- raster::raster(paste0(getwd(), "/data/GIS/habitat_surface.tif"))



# plot original NLCD
plot_nlcd <- rasterVis::levelplot(nlcd, main = "NLCD", att = "class", 
                                  col.regions = nlcd@data@attributes[[1]]$color) + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="black", size = 2))

# plot reclassified NLCD
plot_nlcd_reclassified <- rasterVis::levelplot(nlcd_reclassified, att = "class", 
                                               col.regions = c("darkgreen", "greenyellow", "gray"),
                                               main = "NLCD reclassified") + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="black", size = 2))

# plot original vs reclassified NLCD
plot_nlcd_overall <- gridExtra::grid.arrange(plot_nlcd, 
                                             plot_nlcd_reclassified, ncol = 2)

# plot habitat surface
plot_surface <- rasterVis::levelplot(habitat_surface, main = "Surface metrics", 
                                     margin = FALSE) + 
  latticeExtra::layer(sp::sp.points(sampling_points, pch = 19, col ="black", size = 2))

# plot habitat surface vs. reclassified NLCD``
plot_surface_patch <- gridExtra::grid.arrange(plot_nlcd_reclassified, 
                                              plot_surface, ncol = 2)
