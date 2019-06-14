#### load libraries ####
library(gdistance)
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(lme4)
library(MuMIn)
library(raster)
library(spdep)
library(tidyverse)

# source modular functions
source("scripts/0_modular_functions.R")

#### load data ####

# habitat surface
habitat_surface <- raster::raster("data/GIS/Habitat_surface.tif")

# Rst 
rst <- readr::read_rds("data/rst.rds")

# sample locations
sites <- raster::shapefile("data/GIS/SSR_17_sites.shp")

# conductance values
conductance_surface <- raster::raster("data/GIS/tv_cond.tif")

#### Pre-processing of data ####

# create ids
site_ids <- helpeR::expand_grid_unique(x = seq_along(sites), 
                                       y = seq_along(sites)) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>% 
  purrr::set_names("site_1", "site_2")

# Create ZZ matrix
ZZ_matrix <- create_ZZ(data = site_ids, col_names = c("site_1", "site_2"))

# make sure no negative values are present by addind minimum value
habitat_surface[] <- habitat_surface[] + ceiling(abs(min(habitat_surface[], 
                                                         na.rm = TRUE)))

# calculate transition path
transition_conductance_surface <- gdistance::transition(x = conductance_surface, 
                                                        transitionFunction = mean, directions = 8)

# geographic correction
transition_conductance_surface <- gdistance::geoCorrection(x = transition_conductance_surface,
                                                           type = "c", multpl = FALSE)

#### plot least cost paths ####
# # get neighbouring sampling points
# neighbours <- spdep::tri2nb(sites@coords, row.names = sites$SiteName)
# 
# # plot paths between neighbours
# plot(raster::raster(transition_cost_surface), 
#      xlab = "x coordinate (m)", 
#      ylab = "y coordinate (m)", 
#      legend.lab = "Conductance")
# 
# # plot sampling points
# points(sites)
# 
# # plot neighbouring lines
# plot(neighbours, sites@coords, col = "darkgrey", add = TRUE)
# 
# # add paths
# for (focal_id in seq_along(neighbours)) {
#   
#   # get all neighbours
#   neighbours_focal <- neighbours[[focal_id]]
#   
#   # remove neighbours already considered
#   neighbours_focal <- neighbours_focal[neighbours_focal > focal_id]
#   
#   for (other_id in seq_along(neighbours_focal)) {
#     
#     # get neighbouring id
#     neighbours_other <- neighbours_focal[[other_id]]
#     
#     # calculate path
#     shortest_path <- gdistance::shortestPath(x = transition_cost_surface, 
#                                              origin = sites[focal_id, ], 
#                                              goal = sites[neighbours_other, ], 
#                                              output = "SpatialLines")
#     
#     # plot line
#     lines(shortest_path, col = "green", lwd = 1.5)
#   }
# }

# calculate least-cost distances
distance_least_cost <- gdistance::costDistance(x = transition_conductance_surface,
                                               fromCoords = sites)

# save results
helpeR::save_rds(object = distance_least_cost, 
                 filename = "distance_least_cost.rds", 
                 path = "data/Output/", overwrite = FALSE)

# read already computed data
# distance_least_cost <- readr::read_rds("data/Output/distance_least_cost.rds")

distance_least_cost_df <- tibble::tibble(site_1 = site_ids$site_1, 
                                         site_2 = site_ids$site_2, 
                                         least_cost = as.numeric(distance_least_cost)) %>% 
  dplyr::left_join(rst, by = c("site_1", "site_2"))

# calculate resistance distances (this takes some time)
distance_resistance <- gdistance::commuteDistance(x = transition_conductance_surface,
                                                  coords = sites)

# save results
helpeR::save_rds(object = distance_resistance, 
                 filename = "distance_resistance.rds", 
                 path = "data/Output/", overwrite = FALSE)

# read already computed data
# distance_resistance <- readr::read_rds("data/Output/distance_resistance.rds")

# save distances in tibble
distance_resistance_df <- tibble::tibble(site_1 = site_ids$site_1, 
                                         site_2 = site_ids$site_2,
                                         resistance = as.numeric(distance_resistance)) %>% 
  dplyr::left_join(rst, by = c("site_1", "site_2"))

# calculate correlation between distances
distance_correlation <- cor(x = distance_least_cost_df$least_cost,
                            y = distance_resistance_df$resistance,
                            method = "spearman")

plot(distance_least_cost_df$least_cost ~ distance_resistance_df$resistance)

#### Run models #### 

# least distance cost model
# MH: I think we want REML = TRUE since we are not dredging the model
modell_least_cost_REML <- modular_function(variables = RST ~ least_cost + (1|site_1), 
                                           data = distance_least_cost_df, 
                                           REML = TRUE, 
                                           ZZ = ZZ_matrix)

summary(modell_least_cost_REML)
MuMIn::r.squaredGLMM(modell_least_cost_REML)
MuMIn::AICc(modell_least_cost_REML)

# resistance model
modell_resistance_REML <- modular_function(variables = RST ~ resistance + (1|site_1),
                                           data = distance_resistance_df,
                                           REML = TRUE, 
                                           ZZ = ZZ_matrix)

summary(modell_resistance_REML)
MuMIn::r.squaredGLMM(modell_resistance_REML)
MuMIn::AICc(modell_resistance_REML)
