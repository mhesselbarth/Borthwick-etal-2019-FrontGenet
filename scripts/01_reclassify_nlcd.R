library(gridExtra)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt") - only needed to save data
library(raster)
library(rasterVis)
library(tidyverse)

# import GIS layer
nlcd <- raster::raster(paste0(getwd(), "/data/GIS/nlcd_in.tif"))

# read raster to memory
nlcd <- raster::readAll(nlcd)

# convert values to factor
nlcd <- raster::ratify(nlcd)

# show all unique landcover classes 
# sort(unique(raster::values(nlcd)))

# add classes of NLCD
levels_nlcd <- raster::levels(nlcd)[[1]]
levels_nlcd$class <- c("water", "urban", "residential (lo)", "residential (hi)", "commercial", 
                       "rock", "forest (d)", "forest (e)", "forest (m)", "shrub", "grass",
                       "pasture", "crops", "wetlands (w)", "wetlands (h)", "NA")

# add RAT to raster
levels(nlcd) <- levels_nlcd

# new raster to reclassify
nlcd_reclassified <- nlcd

# reclassify land-cover classes
nlcd_reclassified[nlcd_reclassified %in% c(41, 42, 43)] <- 1 # forest habitat
nlcd_reclassified[nlcd_reclassified %in% c(52, 71, 81)] <- 2 # complementary habitat
nlcd_reclassified[!nlcd_reclassified %in% c(1, 2) & 
                    !is.na(nlcd_reclassified)] <- 3 # non-habitat

# convert values to factor
nlcd_reclassified <- raster::ratify(nlcd_reclassified)

# add classes of reclassification
levels_nlcd_reclassified <- levels(nlcd_reclassified)[[1]]
levels_nlcd_reclassified$class <- c("forest", "complementary", "non-habitat")

# add RAT to raster
levels(nlcd_reclassified) <- levels_nlcd_reclassified

# Save results
suppoRt::save_rds(object = nlcd, 
                  filename = "nlcd.rds", 
                  path = paste0(getwd(), "/data/output"), 
                  overwrite = FALSE)

# Save results
suppoRt::save_rds(object = nlcd_reclassified, 
                  filename = "nlcd_reclassified.rds", 
                  path = paste0(getwd(), "/data/output"), 
                  overwrite = FALSE)
