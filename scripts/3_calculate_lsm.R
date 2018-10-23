# load libraries
library(dplyr)
library(landscapemetrics)
library(magrittr)
library(purrr)
library(raster)
library(rgdal)
library(stringr)

# load the clippings
clippings_pmm <- readRDS(paste0(getwd(), "/data/output/clippings_pmm.rds"))

names_clippings <- purrr::map_chr(clippings_pmm, function(x) names(x))

names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE)

# Calculate landscape-level metrics
# metrics <- c("lsm_l_lpi", "lsm_l_area_mn", "lsm_l_split")

landscape_metrics <- landscapemetrics::calculate_lsm(clippings_pmm, 
                                                     what = "landscape")

# Add name of sites
landscape_metrics <- dplyr::mutate(landscape_metrics, 
                                   site_a = as.integer(names_clippings[layer, 2]), 
                                   site_b = as.integer(names_clippings[layer, 3]))

# Import sample points
sample_points <- getwd() %>% 
  paste0("/data/SSR_17_Sites.shp") %>%
  rgdal::readOGR() %>%
  tibble::as.tibble()

# Calculate distance between each point
distance_matrix <- sample_points %>%
  dplyr::select(E, N) %>% # Only select the coordinates
  dist(diag = TRUE, upper = TRUE) %>%
  as.matrix()

# Add euclidean distance to each pair of sitres
landscape_metrics <- dplyr::mutate(landscape_metrics, euclidean_distance = distance_matrix[cbind(site_a, site_b)]) %>% 
  dplyr::arrange(site_a, site_b)

# Order and save results
write.table(landscape_metrics, # Write to comma separeted file
            file = paste0(getwd(), '/data/output/landscape_metrics.csv'), 
            sep = ";", dec = ".",
            row.names = FALSE)
                 