# library(ggsn)
# library(patchwork)
library(raster)
library(sf)
library(sp)
library(tidyverse)

#### General data ####

# read state boundary
boundary <- sf::read_sf("data/TV_maps_data/IN_poly.shp")

# convert do data frame for easier plotting
boundary_df <- sf::st_coordinates(boundary) %>%
  tibble::as_tibble()

# read sample sites
sampling_sites <- sf::read_sf("data/TV_maps_data/TV_sites.shp")

# plot sample sites and state border
sampling_sites_gg <- ggplot2::ggplot(data = sampling_sites) + 
  ggplot2::geom_point(ggplot2::aes(x = E, y = N, col = Site_name), size = 2, pch = 19) + 
  ggplot2::geom_polygon(data = boundary_df, 
                        ggplot2::aes(x = X, y = Y), fill = NA, col = "black") + 
  ggplot2::scale_color_viridis_d(name = "Sampling sites") +
  ggplot2::coord_equal() + 
  ggplot2::theme_void()

# save plot
ggplot2::ggsave("plots/sampling_sites_gg.png", plot = sampling_sites_gg)

#### Habitat suitability surface #### 

# import data
habitat_surface <- raster::raster("data/TV_maps_data/Habitat_surface.tif")

# convert to data frame
habitat_surface_df <- as.data.frame(habitat_surface, xy = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(Habitat_surface))

# plot habitat surface
habitat_surface_gg <- ggplot2::ggplot(data = habitat_surface_df) + 
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = Habitat_surface)) + 
  ggplot2::geom_polygon(data = boundary_df, 
                        ggplot2::aes(x = X, y = Y), fill = NA, col = "black") + 
  ggplot2::scale_fill_viridis_c("Habitat suitability surface") + 
  ggplot2::coord_equal() + 
  ggplot2::theme_void()

# make sure RAM gets empty
rm(list = c("habitat_surface", "habitat_surface_df"))
gc()

# save plot
ggplot2::ggsave("plots/habitat_surface_gg.png", plot = habitat_surface_gg)

#### Conductance/Resistance surface #### 

# import data 
resistance_surface <- raster::raster("data/TV_maps_data/TV_resist.tif")

# convert to data frame
resistance_surface_df <- as.data.frame(resistance_surface, xy = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(TV_resist))

# plot resistance surface
resistance_surface_gg <- ggplot2::ggplot(data = resistance_surface_df) + 
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = TV_resist)) + 
  ggplot2::geom_polygon(data = boundary_df, 
                        ggplot2::aes(x = X, y = Y), fill = NA, col = "black") + 
  ggplot2::scale_fill_viridis_c("Resistance surface") + 
  ggplot2::coord_equal() + 
  ggplot2::theme_void()

# make sure RAM gets empty
rm(list = c("resistance_surface", "resistance_surface_df"))
gc()

# save plot
ggplot2::ggsave("plots/resistance_surface_gg.png", plot = resistance_surface_gg)

#### Current density map #### 

# import data 
current_density <- raster::raster("data/TV_maps_data/TV_current_dens.tif")

# mask data to Indiana
current_density <- raster::mask(x = current_density, mask = as(boundary, "Spatial"))

# convert to data frame
current_density_df <- as.data.frame(current_density, xy = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(TV_current_dens))

current_density_gg <- ggplot2::ggplot(data = current_density_df) + 
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = TV_current_dens)) + 
  ggplot2::geom_polygon(data = boundary_df, 
                        ggplot2::aes(x = X, y = Y), fill = NA, col = "black") + 
  ggplot2::scale_fill_viridis_c("Current density") + 
  ggplot2::coord_equal() + 
  ggplot2::theme_void()

# make sure RAM gets empty
rm(list = c("current_density", "current_density_df"))
gc()

# save plot
ggplot2::ggsave("plots/current_density_gg.png", plot = current_density_gg)

#### Reclassified habitat patches map ####

# import data 
nlcd <- readr::read_rds("data/Output/nlcd_reclassified.rds")

# convert to data frame (this takes quite some time)
nlcd_df <- as.data.frame(nlcd, xy = TRUE) %>%
  tibble::as_tibble() %>%
  dplyr::filter(!is.na(layer_class))

# plot NLCD
nlcd_gg <- ggplot2::ggplot(data = nlcd_df) + 
  ggplot2::geom_raster(ggplot2::aes(x = x, y = y, fill = layer_class)) + 
  ggplot2::geom_polygon(data = boundary_df, 
                        ggplot2::aes(x = X, y = Y), fill = NA, col = "black") + 
  
  ggplot2::scale_fill_manual(name = "reclassified NLCD", 
                             values = c("forest" = "darkgreen", 
                                        "complementary" = "greenyellow", 
                                        "non-habitat" = "gray")) +
  ggplot2::coord_equal() + 
  ggplot2::theme_void()

# make sure RAM gets empty
rm(list = c("nlcd", "nlcd_df"))
gc()

# save plots
ggplot2::ggsave("plots/nlcd_gg.png", plot = nlcd_gg)

#### Combine plots ####
# maps_overall_gg <- habitat_surface_gg + resistance_surface_gg + current_density_gg + nlcd_gg + plot_layout(ncol = 2, nrow = 2)
# 
# ggplot2::ggsave("plots/maps_overall_gg.png", plot = maps_overall_gg)
