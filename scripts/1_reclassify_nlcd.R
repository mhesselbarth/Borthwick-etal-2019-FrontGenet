library(gridExtra)
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions") - only needed to save data
library(raster)
library(rasterVis)
library(tidyverse)

# Import GIS layer
nlcd <- raster::raster(paste0(getwd(), "/data/GIS/nlcd_in.tif"))

# show all unique landcover classes 
# sort(unique(raster::values(nlcd)))

# new raster to reclassify
nlcd_reclassified <- nlcd

# reclassify land-cover classes
# forest habitat
raster::values(nlcd_reclassified)[raster::values(nlcd_reclassified) %in% c(41, 42, 43)] <- 1
# complementary habitat
raster::values(nlcd_reclassified)[raster::values(nlcd_reclassified) %in% c(52, 71, 81)] <- 2
# non-habitat
raster::values(nlcd_reclassified)[!raster::values(nlcd_reclassified) %in% c(1, 2) & 
                                    !is.na(raster::values(nlcd_reclassified))] <- 3

# convert to factor
raster::values(nlcd) <- factor(raster::values(nlcd), 
                               levels = c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95, 127), 
                               labels = c("water", "urban", "residential (lo)", "residential (hi)", "commercial", 
                                          "rock", "forest (d)", "forest (e)", "forest (m)", "shrub", "grass",
                                          "pasture", "crops", "wetlands (w)", "wetlands (h)", "NA"))

# convert to factor
raster::values(nlcd_reclassified) <- factor(raster::values(nlcd_reclassified), 
                                            levels = c(1, 2, 3), 
                                            labels = c("forest", "complementary", "non-habitat"))

# Save results
UtilityFunctions::save_rds(object = nlcd, 
                           filename = "nlcd.rds", 
                           path = paste0(getwd(), "/data/output"), 
                           overwrite = FALSE)

# Save results
UtilityFunctions::save_rds(object = nlcd_reclassified, 
                           filename = "nlcd_reclassified.rds", 
                           path = paste0(getwd(), "/data/output"), 
                           overwrite = FALSE)
