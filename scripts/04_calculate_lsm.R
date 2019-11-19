# load libraries
library(clustermq)
library(suppoRt) # devtools::install_github("mhesselbarth/helpeR")
library(landscapemetrics)
library(raster)
library(sp)
library(tidyverse)

source(paste0(getwd(), "/scripts/00_calculate_lsm_helper.R"))
source(paste0(getwd(), "/scripts/00_clip_and_calc.R"))

#### Load data ####

# # load the clippings
# clippings_pmm_nlcd <- readRDS(paste0(getwd(), "/data/Output/clippings_pmm_nlcd.rds"))
# 
# # check if all rasters all loaded in memory
# all(purrr::map_lgl(clippings_pmm_nlcd, raster::inMemory))
# 
# # extract names
# names_clippings <- purrr::map_chr(clippings_pmm_nlcd, function(x) names(x))
# 
# names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

# # load input layer
# nlcd_layer <- readRDS(paste0(getwd(), "/data/Output/nlcd_reclassified.rds"))

# load sampling points
sampling_points <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))

#### Specify metrics ####
landscape_sub <- c("lsm_l_ai", 
                   "lsm_l_area_mn", 
                   "lsm_l_cai_mn", 
                   "lsm_l_condent", 
                   "lsm_l_contag", 
                   "lsm_l_core_mn", 
                   "lsm_l_division", 
                   "lsm_l_ed", 
                   "lsm_l_ent", 
                   "lsm_l_iji", 
                   "lsm_l_joinent",
                   "lsm_l_lpi", 
                   "lsm_l_lsi", 
                   "lsm_l_mesh", 
                   "lsm_l_mutinf", 
                   "lsm_l_np", 
                   "lsm_l_pd", 
                   "lsm_l_pladj", 
                   "lsm_l_pr", 
                   "lsm_l_prd", 
                   "lsm_l_rpr", 
                   "lsm_l_shdi", 
                   "lsm_l_shei",
                   "lsm_l_sidi",
                   "lsm_l_siei",
                   "lsm_l_split", 
                   "lsm_l_ta", 
                   "lsm_l_te")

# #### Calculate locally ####
# # Calculate metrics locally
# landscape_metrics <- landscapemetrics::calculate_lsm(clippings_pmm,
#                                                      what = landscape_sub,
#                                                      classes_max = 3)
# 
# # Calculate metrics locally but overall printing progress
# total_clippigings <- length(clippings_pmm_nlcd)
# 
# landscape_metrics <- purrr::map(seq_along(clippings_pmm_nlcd), function(x) {
# 
#   print(paste0("Progress: ", x, " from ", total_clippigings))
# 
#   result <- calculate_lsm(landscape = clippings_pmm_nlcd[[x]],
#                           what = landscape_sub,
#                           classes_max = 3,
#                           verbose = FALSE,
#                           progress = FALSE)
# 
#   gc(verbose = FALSE, reset = TRUE, full = TRUE)
# 
#   return(result)
# })

#### clustermq (clip_and_calc) ####

# get all combinations
sampling_ids <- suppoRt::expand_grid_unique(x = seq_along(sampling_points),
                                            y = seq_along(sampling_points))

# run metrics
landscape_metrics <- suppoRt::submit_to_cluster(fun = clip_and_calc,
                                                focal_plot = sampling_ids[, 1],
                                                other_plot = sampling_ids[, 2],
                                                n_jobs = nrow(sampling_ids),
                                                log_worker = TRUE,
                                                const = list(sampling_points = sampling_points,
                                                             # input_layer = nlcd_layer,
                                                             what = landscape_sub,
                                                             classes_max = 3, 
                                                             path = "/home/uni08/hesselbarth3/nlcd_reclassified.rds"),
                                                template = list(queue = "medium",
                                                                walltime = "02:00:00",
                                                                mem_cpu = "8192",
                                                                processes = 1))

suppoRt::save_rds(object = landscape_metrics,
                  filename = "landscape_metrics_raw.rds",
                  path = paste0(getwd(), "/data/Output"),
                  overwrite = FALSE)

# bind to one dataframe 
landscape_metrics <- dplyr::bind_rows(landscape_metrics)

# replace layer with 1:136
landscape_metrics$layer <- rep(x = 1:nrow(sampling_ids), 
                               each = length(unique(landscape_metrics$metric)))

suppoRt::save_rds(object = landscape_metrics,
                  filename = "landscape_metrics.rds",
                  path = paste0(getwd(), "/data/Output"),
                  overwrite = FALSE)
