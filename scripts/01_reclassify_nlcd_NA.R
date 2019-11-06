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
nlcd_reclassified_NA <- nlcd

# reclassify land-cover classes
nlcd_reclassified_NA[nlcd_reclassified_NA %in% c(41, 42, 43)] <- 1 # forest habitat
nlcd_reclassified_NA[nlcd_reclassified_NA %in% c(52, 71, 81)] <- 2 # complementary habitat
nlcd_reclassified_NA[!nlcd_reclassified_NA %in% c(1, 2) & 
                    !is.na(nlcd_reclassified_NA)] <- NA # non-habitat

# convert values to factor
nlcd_reclassified_NA <- raster::ratify(nlcd_reclassified_NA) # NA dataset

# add classes of reclassification NA datset
levels_nlcd_reclassified_NA <- levels(nlcd_reclassified_NA)[[1]]
levels_nlcd_reclassified_NA$class <- c("forest", "complementary")

# add RAT to raster
levels(nlcd_reclassified_NA) <- levels_nlcd_reclassified_NA # NA dataset

# Save results
suppoRt::save_rds(object = nlcd, 
                  filename = "nlcd.rds", 
                  path = paste0(getwd(), "/data/output"), 
                  overwrite = FALSE)

# Save results NA dataset
suppoRt::save_rds(object = nlcd_reclassified_NA, 
                  filename = "nlcd_reclassified_NA.rds", 
                  path = paste0(getwd(), "/data/output"), 
                  overwrite = FALSE)
