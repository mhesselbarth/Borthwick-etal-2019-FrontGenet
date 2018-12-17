# load libraries
library(clustermq)
library(landscapemetrics)
library(raster)
library(tidyverse)

source(paste0(getwd(), "/scripts/0_calculate_lsm_helper.R"))

# load the clippings
clippings_pmm <- readRDS(paste0(getwd(), "/data/output/clippings_pmm_nlcd.rds"))

names_clippings <- purrr::map_chr(clippings_pmm, function(x) names(x))
names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

metrics <- landscapemetrics::list_lsm(level = c("class", "landscape"), simplify = TRUE)

metrics <- metrics[!metrics %in% c("lsm_c_contig_mn", 
                                   "lsm_c_contig_sd", 
                                   "lsm_c_contig_cv",
                                   "lsm_l_contig_mn", 
                                   "lsm_l_contig_sd", 
                                   "lsm_l_contig_cv")]

# # Calculate landscape-level metrics locally
# # Nicer code but can't print overall progress at the moment...
# landscape_metrics <- landscapemetrics::calculate_lsm(clippings_pmm,
#                                                      what = c("class", "landscape"))
# 
# # Add name of sites
# landscape_metrics <- dplyr::mutate(landscape_metrics,
#                                    site_a = as.integer(names_clippings[layer, 2]),
#                                    site_b = as.integer(names_clippings[layer, 3]))
# 
# # Print progress
# landscape_metrics <- purrr::map(seq_along(clippings_pmm), function(x) {
# 
#   cat(paste0("\r> Progress: ", x, " from ", length(clippings_pmm)))
# 
#   landscapemetrics::calculate_lsm(clippings_pmm[[x]], what = c("class", "landscape"),
#                                   verbose = FALSE)
# })
# 
# # Add name of sites
# landscape_metrics <- dplyr::bind_rows(landscape_metrics, .id = "layer_bind_rows") %>%
#   dplyr::mutate(layer = as.integer(layer_bind_rows),
#                 site_a = as.integer(names_clippings[layer, 2]),
#                 site_b = as.integer(names_clippings[layer, 3])) %>%
#   dplyr::select(-layer_bind_rows)
# 
# Calculate landscape-level metrics on high performance cluster
landscape_metrics <- clustermq::Q(fun = calculate_lsm_helper,
                                  landscape = clippings_pmm,
                                  const = list(what = metrics,
                                               classes_max = 3),
                                  n_jobs = length(clippings_pmm),
                                  template = list(queue = "mpi-long",
                                                  walltime = "120:00",
                                                  processes = 1))

# Rowbind returning list and add site names
landscape_metrics <- dplyr::bind_rows(landscape_metrics, .id = "layer_bind_rows") %>%
  dplyr::mutate(layer = as.integer(layer_bind_rows),
                site_a = as.integer(names_clippings[layer, 2]),
                site_b = as.integer(names_clippings[layer, 3])) %>%
  dplyr::select(-layer_bind_rows)

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
UtilityFunctions::save_rds(object = landscape_metrics, 
                          filename = "landscape_metrics.rds", 
                          path = paste0(getwd(), "/data/output"), 
                          overwrite = FALSE)

write.table(landscape_metrics,
            file = paste0(getwd(), '/data/output/landscape_metrics.csv'),
            sep = ";", dec = ".",
            row.names = FALSE)
