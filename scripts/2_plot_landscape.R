
# Import libraries
library(patchwork) # devtools::install_github("thomasp85/patchwork") - arrange ggplots next to each other
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions") - only needed to save data
library(raster)
library(sp)
library(tidyverse)

# Import sampling points
sampling_points <- getwd() %>% # path to folder containing GIS data
  paste0("/data/GIS/SSR_17_sites.shp") %>% # string containing path and name of shape file
  raster::shapefile() # read data
names(sampling_points) <- c("x", "y", "Group", "Site_name") # rename columns

 # Import state boarder Indiana
boarder_indiana <- getwd() %>% # path to folder containing GIS data
  paste0("/data/GIS/indiana_state.shp") %>% # string containing path and name of shape file
  raster::shapefile() %>% # read data
  sp::spTransform("+proj=utm +zone=16 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") %>% # set CRS
  raster::aggregate(dissolve = T) # polygon with only the state boarder of Indiana 

# Import GIS layer
# forest_300 <- getwd() %>% 
#   paste0("/data/GIS/IN_forest_300m.tif") %>%
#   raster::raster()
# 
# forest_300_relative <- getwd() %>% 
#     paste0("/data/GIS/IN_forest_2.1Wind_300m.tif") %>%
#     raster::raster()
# 
# forest_cover_30 <- getwd() %>% 
#     paste0("/data/GIS/IN_forest_cover_30m.tif") %>%
#     raster::raster()
# 
# # Set NA to 0 
# forest_300[is.na(forest_300)] <- 0 # Set NA values to 0, because there are NAs inside the state
# forest_300 <- raster::mask(x = forest_300, mask = boarder_indiana) # clip/mask data only to the state boarder

# Import GIS layer
nlcd <- getwd() %>%
  paste0("/data/GIS/nlcd_in.tif") %>%
  raster::raster()

# convert raster to dataframe with xy-coordinates and z-values
nlcd_df <- raster::as.data.frame(nlcd_df, 
                                 na.rm = TRUE, xy = TRUE)
  

# Import habitat surface
habitat_surface <- getwd() %>% 
  paste0("/data/GIS/habitat_surface.tif") %>%
  raster::raster()

# convert raster to dataframe with xy-coordinates and z-values
habitat_surface_df <- raster::as.data.frame(habitat_surface, 
                                            na.rm = TRUE, xy = TRUE)
# Plot gis layer with sampling points
ggplot_nlcd_df <- ggplot2::ggplot(data = nlcd_df) + # plot data
    ggplot2::geom_raster(aes(x = x, y = y, fill = factor(IN_forest_300m))) + # plot patches with fill according to value
    ggplot2::geom_point(data = tibble::as.tibble(sampling_points), # plot sampling points
                        aes(x = x, y = y, col = "red"), shape = 19, size = 2) +
    ggplot2::coord_equal() + # ratio of 1:1 of x and y axis
    ggplot2::scale_fill_viridis_d(name = "Landcover", labels = c("Non-forest", "Forest")) + # setting colours of fill
    ggplot2::scale_color_manual(name = "", values = "red", labels = "Sampling points") + # set labels legend
    ggplot2::labs(x = "x coordinate", y = "y-coordinate", title = "Forest cover - 300m") + # set labels axis
    ggplot2::theme_classic() + # overall appearance of the plot
    ggplot2::theme(axis.text=element_blank(), # remove xy-coordinates from axis (not very informative)
                   axis.ticks=element_blank())

ggplot_habitat_surface <- ggplot2::ggplot(data = habitat_surface_df) + # plot data
  ggplot2::geom_raster(aes(x = x, y = y, fill = habitat_surface)) + # plot patches with fill according to value
  ggplot2::geom_point(data = tibble::as.tibble(sampling_points), # plot sampling points
                      aes(x = x, y = y, col = "red"), shape = 19, size = 2) +
  ggplot2::coord_equal() + # ratio of 1:1 of x and y axis
  ggplot2::scale_fill_viridis_c(name = "Habitat quality") + # setting colours of fill
  ggplot2::scale_color_manual(name = "", values = "red", labels = "Sampling points") + # set labels legend
  ggplot2::labs(x = "x coordinate", y = "y-coordinate", title = "Habitat quality") + # set labels axis
  ggplot2::theme_classic() + # overall appearance of the plot
  ggplot2::theme(axis.text=element_blank(), # remove xy-coordinates from axis (not very informative)
                 axis.ticks=element_blank())

ggplot_surface_forest <- ggplot_habitat_surface + ggplot_forest

# Save results
overwrite <- FALSE

UtilityFunctions::save_ggplot(plot = ggplot_surface_forest, 
                              filename = "ggplot_surface_forest.png", 
                              path = paste0(getwd(), "/plots"), 
                              overwrite = overwrite)

