# load libraries
# library(clustermq)
library(landscapemetrics)
library(raster)
library(tidyverse)

source(paste0(getwd(), "/scripts/0_calculate_lsm_helper.R"))

# load the clippings
clippings_pmm <- readRDS(paste0(getwd(), "/data/output/clippings_pmm.rds"))

names_clippings <- purrr::map_chr(clippings_pmm, function(x) names(x))
names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

# Calculate landscape-level metrics locally
landscape_metrics <- landscapemetrics::calculate_lsm(clippings_pmm,
                                                     what = "landscape")

# Add name of sites
landscape_metrics <- dplyr::mutate(landscape_metrics,
                                   site_a = as.integer(names_clippings[layer, 2]),
                                   site_b = as.integer(names_clippings[layer, 3]))

# Calculate landscape-level metrics on high performance cluster
# landscape_metrics <- clustermq::Q(fun = calculate_lsm_helper,
#                                   landscape = clippings_pmm, 
#                                   const = list(what = "landscape"), 
#                                   n_jobs = length(clippings_pmm),
#                                   template = list(queue = "mpi", 
#                                                   walltime = "48:00", 
#                                                   processes = 1))
# 
# Rowbind returning list and add site names
# landscape_metrics <- dplyr::bind_rows(landscape_metrics, .id = "layer_bind_rows") %>%
#   dplyr::mutate(layer = as.integer(layer_bind_rows),
#                 site_a = as.integer(names_clippings[layer, 2]),
#                 site_b = as.integer(names_clippings[layer, 3])) %>%
#   dplyr::select(-layer_bind_rows)

# Import sample points
sample_points <- getwd() %>%
  paste0("/data/GIS/SSR_17_sites.shp") %>%
  raster::shapefile() %>%
  tibble::as_tibble()

# Calculate distance between each point
distance_matrix <- sample_points %>%
  dplyr::select(E, N) %>%
  dist(diag = TRUE, upper = TRUE) %>%
  as.matrix()

# Add euclidean distance to each pair of sitres
landscape_metrics <- dplyr::mutate(landscape_metrics,
                                   euclidean_distance = distance_matrix[cbind(site_a, site_b)]) %>% 
  dplyr::arrange(site_a, site_b)

# Order and save results
write.table(landscape_metrics_df,
            file = paste0(getwd(), '/data/output/landscape_metrics.csv'), 
            sep = ";", dec = ".",
            row.names = FALSE)
