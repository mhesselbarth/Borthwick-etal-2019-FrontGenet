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
# distance_lcp <- gdistance::costDistance(x = transition_conductance_surface,
#                                         fromCoords = sites)
# 
# # save results
# helpeR::save_rds(object = distance_lcp, 
#                  filename = "distance_least_cost.rds", 
#                  path = "data/Output/", overwrite = FALSE)

# read already computed data
distance_lcp <- readr::read_rds("data/Output/distance_least_cost.rds")

df_lcp <- tibble::tibble(site_1 = site_ids$site_1, 
                         site_2 = site_ids$site_2, 
                         least_cost = as.numeric(distance_lcp)) %>% 
  dplyr::mutate(least_cost_scaled = as.numeric(scale(least_cost))) %>%
  dplyr::left_join(rst, by = c("site_1", "site_2"))

# calculate resistance distances (this takes some time)
# distance_res <- gdistance::commuteDistance(x = transition_conductance_surface,
#                                            coords = sites)
# 
# # save results
# helpeR::save_rds(object = distance_res, 
#                  filename = "distance_resistance.rds", 
#                  path = "data/Output/", overwrite = FALSE)

# read already computed data
distance_res <- readr::read_rds("data/Output/distance_resistance.rds")

# save distances in tibble
df_res <- tibble::tibble(site_1 = site_ids$site_1,
                         site_2 = site_ids$site_2,
                         resistance = as.numeric(distance_res)) %>% 
  dplyr::mutate(resistance_scaled = as.numeric(scale(resistance))) %>%
  dplyr::left_join(rst, by = c("site_1", "site_2"))

# # calculate correlation between distances
# distance_correlation <- cor(x = df_lcp$least_cost,
#                             y = df_res$resistance,
#                             method = "spearman")
# 
# plot(df_lcp$least_cost ~ df_res$resistance)

#### Run models #### 

# least distance cost model
model_dredge_lcp <- modular_function(response = "RST", 
                                     explanatory = "least_cost", 
                                     random = "(1|site_1)",
                                     data = df_lcp,
                                     ZZ = ZZ_matrix, 
                                     REML = FALSE) %>% 
  get_model_aic(n = 136) %>%
  dplyr::bind_cols(model = 1, .)

r2_lcp <- modular_function(response = "RST", 
                           explanatory = "least_cost",
                           random = "(1|site_1)",
                           data = df_lcp,
                           ZZ = ZZ_matrix, 
                           REML = TRUE) %>% 
  get_model_r2() %>%
  dplyr::bind_cols(model = 1, .)

# resistance model
model_dredge_res <- modular_function(response = "RST", 
                                     explanatory = "resistance_scaled", 
                                     random = "(1|site_1)",
                                     data = df_res,
                                     ZZ = ZZ_matrix,
                                     REML = FALSE) %>%
  get_model_aic(n = 136) %>%
  dplyr::bind_cols(model = 1, .)

r2_res <- modular_function(response = "RST", 
                           explanatory = "resistance_scaled", 
                           random = "(1|site_1)",
                           data = df_res,
                           ZZ = ZZ_matrix,
                           REML = TRUE) %>%
  get_model_r2() %>%
  dplyr::bind_cols(model = 1, .)

#### Compare all models ####
model_dredge_lcp
model_dredge_res

r2_lcp
r2_res
