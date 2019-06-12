#### load libraries ####
library(gdistance)
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(lme4)
library(MuMIn)
library(raster)
library(spdep)
library(tidyverse)

#### load modular functions ####

# load modelling functions
# model optimization function
modular_function <- function(variables, data, REML = TRUE,
                             ZZ = NULL, call = NULL) {
  
  # parse the data and formula
  model <- lme4::lFormula(variables, data = data, REML = REML)
  
  if (!is.null(ZZ)) {
    
    # replace ZZ matrix
    model$reTrms$Zt <- ZZ
    
    message("> Replaced 'model$reTrms$Zt <- ZZ'")
  }
  
  # create the deviance function to be optimized
  deviance_function <- do.call(lme4::mkLmerDevfun, model)
  
  # optimize the deviance function:
  optimization <- lme4::optimizeLmer(deviance_function)
  
  # package up the results
  model <- lme4::mkMerMod(rho = environment(deviance_function),
                          opt = optimization,
                          reTrms = model$reTrms,
                          fr = model$fr)
  
  if (!is.null(call)) {
    
    # replace call for dredging
    model@call <- call
    
    message("> Replaced 'model@call <- call'")
    
  }
  
  return(model)
}

get_model_info <- function(model, n) {
  
  # getting AIC and r2
  # should we use extractAIC here ?
  information_criterion <- purrr::map_dfr(model, function(x) {
    
    ic <- tibble::tibble(AIC = AIC(x),
                         BIC = BIC(x),
                         log_like = attr(logLik(x), "df"))
    
    r2 <- MuMIn::r.squaredGLMM(x)
    
    ic <- dplyr::mutate(ic, 
                        r2_marginal = r2[1, 1], 
                        r2_conditional = r2[1, 2])
  }, .id = "model")
  
  # calculating corrected AICc and BICc and weights
  information_criterion <- dplyr::mutate(information_criterion,
                                         AICc = AIC + 2 * log_like *
                                           (log_like + 1) / (n - log_like - 1),
                                         AICc_ew = exp(-0.5 * (AICc - min(AICc))) /
                                           sum(exp(-0.5 * (AICc - min(AICc)))),
                                         BIC_ew = exp(-0.5 * (BIC - min(BIC))) /
                                           sum(exp(-0.5 * (BIC - min(BIC)))))
  
  return(information_criterion)
}

#### load data ####

# habitat surface
habitat_surface <- raster::raster("data/GIS/habitat_surface.tif")

# Rst 
rst <- readr::read_rds("data/rst.rds")

# sample locations
sites <- raster::shapefile("data/GIS/SSR_17_sites.shp")

# conductance values
conductance_surface <- raster::raster("data/TV_maps_data/tv_cond.tif")


#### Pre-processing of data ####

# create ids
site_ids <- helpeR::expand_grid_unique(x = seq_along(sites), 
                                       y = seq_along(sites))

# make sure no negative values are present by addind minimum value
habitat_surface[] <- habitat_surface[] + ceiling(abs(min(habitat_surface[], 
                                                         na.rm = TRUE)))

# create cost surface
# MH: Why are we never using this again?
# cost_surface <- 1 / habitat_surface 

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

distance_least_cost_df <- tibble::tibble(site_1 = site_ids[, 1], 
                                         site_2 = site_ids[, 2], 
                                         least_cost = as.numeric(distance_least_cost)) %>% 
  dplyr::left_join(rst, by = c("site_1", "site_2"))

# # calculate resistance distances (this takes some time)
distance_resistance <- gdistance::commuteDistance(x = transition_conductance_surface,
                                                  coords = sites)

# save distances in tibble
distance_resistance_df <- tibble::tibble(site_1 = site_ids[, 1], site_2 = site_ids[, 2],
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
                                           REML = TRUE)

summary(modell_least_cost_REML)
MuMIn::r.squaredGLMM(modell_least_cost_REML)
MuMIn::AICc(modell_least_cost_REML)

# # resistance model
modell_resistance_REML <- modular_function(variables = RST ~ resistance + (1|site_1),
                                           data = distance_resistance_df,
                                           REML = TRUE)

summary(modell_resistance_REML)
MuMIn::r.squaredGLMM(modell_resistance_REML)
MuMIn::AICc(modell_resistance_REML)

#### Calculate AIC weights #### 

# x <- c(exp(-0.5*0),exp(-0.5*8.6),exp(-0.5*15.2),exp(-0.5*13.1),exp(-0.5*2.1))
# x[1]/sum(x)
# x[2]/sum(x)
# x[3]/sum(x)
# x[4]/sum(x)
# x[5]/sum(x)

