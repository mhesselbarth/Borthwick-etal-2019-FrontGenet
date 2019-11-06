# load libraries
library(clustermq)
library(suppoRt) # devtools::install_github("mhesselbarth/helpeR")
library(landscapemetrics)
library(raster)
library(sp)
library(tidyverse)

source(paste0(getwd(), "/scripts/0_calculate_lsm_helper.R"))
source(paste0(getwd(), "/scripts/0_clip_and_calc.R"))

#### Load data ####

# # load the clippings
# clippings_pmm <- readRDS(paste0(getwd(), "/data/output/clippings_pmm_nlcd.rds"))

# # check if all rasters all loaded in memory
# all(purrr::map_lgl(clippings_pmm, raster::inMemory))
# 
# # extract names
# names_clippings <- purrr::map_chr(clippings_pmm, function(x) names(x))
# 
# names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

# load input layer
nlcd_layer <- readRDS(paste0(getwd(), "/data/Output/nlcd_reclassified.rds"))

# load sampling points
sampling_points <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))

# specify metrics
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
                   "lsm_l_sidi,",
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
# total_clippigings <- length(clippings_pmm)
# 
# landscape_metrics <- purrr::map(seq_along(clippings_pmm), function(x) {
# 
#   print(paste0("Progress: ", x, " from ", total_clippigings))
# 
#   result <- calculate_lsm(landscape = clippings_pmm[[x]],
#                           what = landscape_sub,
#                           classes_max = 3,
#                           verbose = FALSE,
#                           progress = FALSE)
# 
#   gc(verbose = FALSE, reset = TRUE, full = TRUE)
# 
#   return(result)
# })
# 
# #### clustermq ####
# Calculate metrics on high performance cluster
# landscape_metrics <- suppoRt::submit_to_cluster(fun = calculate_lsm_helper,
#                                                 landscape = clippings_pmm,
#                                                 const = list(what = landscape_sub,
#                                                              classes_max = 3),
#                                                 n_jobs = length(clippings_pmm),
#                                                 template = list(queue = "fat",
#                                                                 walltime = "06:00:00",
#                                                                 processes = 1))
# 
# suppoRt::save_rds(object = landscape_metrics,
#                   filename = "landscape_metrics_raw.rds",
#                   path = paste0(getwd(), "/data/output"),
#                   overwrite = TRUE)
# 
# # Rowbind returning list and add site names
# landscape_metrics <- dplyr::bind_rows(landscape_metrics, .id = "layer_bind_rows") %>%
#   dplyr::mutate(layer = as.integer(layer_bind_rows),
#                 site_a = as.integer(names_clippings[layer, 2]),
#                 site_b = as.integer(names_clippings[layer, 3])) %>%
#   dplyr::select(-layer_bind_rows)
# 
# # Import sample points
# sample_points <- getwd() %>%
#   paste0("/data/GIS/SSR_17_sites.shp") %>%
#   raster::shapefile() %>%
#   tibble::as_tibble()
# 
# # Calculate distance between each point
# distance_matrix <- sample_points %>%
#   dplyr::select(E, N) %>%
#   dist(diag = TRUE, upper = TRUE) %>%
#   as.matrix()
# 
# # Add euclidean distance to each pair of sitres
# landscape_metrics <- dplyr::mutate(landscape_metrics,
#                                    euclidean_distance = distance_matrix[cbind(site_a, 
#                                                                               site_b)]) %>% 
#   dplyr::arrange(site_a, site_b)
# 
# # Order and save results
# suppoRt::save_rds(object = landscape_metrics, 
#                   filename = "landscape_metrics.rds", 
#                   path = paste0(getwd(), "/data/output"), 
#                   overwrite = FALSE)
# 
# write.table(landscape_metrics,
#             file = paste0(getwd(), '/data/output/landscape_metrics.csv'),
#             sep = ";", dec = ".",
#             row.names = FALSE)
# 
# #### FUTURE for HPC ####
# 
# # load the packages
# library("future")
# library("future.batchtools")
# library("furrr")
# 
# # now we specify a future topology that fits our HPC
# # login node -> cluster nodes -> core/ multiple cores
# login <- tweak(future::remote,
#                workers = "gwdu103.gwdg.de",
#                user = "hesselbarth3") # user = login credential
# 
# sbatch <- tweak(future.batchtools::batchtools_slurm,
#                 template = "future_slurm.tmpl",
#                 resources = list(job.name = "calculate_lsm", # name of the job
#                                  log.file = "calculate_lsm.log", # name of log file
#                                  queue = "fat", # which partition
#                                  service = "normal", # which QOS
#                                  walltime = "01:00:00", # walltime <hh:mm:ss>
#                                  processes = 1)) # number of cores
# 
# future::plan(list(login, sbatch, future::sequential)) # how to run on nodes, could also be sequential
# 
# # no max size of globals
# options(future.globals.maxSize = Inf)
# 
# landscape_metrics %<-% furrr::future_map(clippings_pmm, function(x) {
# 
#   calculate_lsm(landscape = x,
#                 what = landscape_sub,
#                 classes_max = 3,
#                 verbose = FALSE,
#                 progress = FALSE)
# })
# 
# future::resolved(future::futureOf(landscape_metrics))
# 
#### clustermq (clip_and_calc) ####

# get all combinations
sampling_ids <- suppoRt::expand_grid_unique(x = seq_along(sampling_points),
                                            y = seq_along(sampling_points))

# run metrics
landscape_metrics <- suppoRt::submit_to_cluster(fun = clip_and_calc,
                                                focal_plot = sampling_ids[, 1],
                                                other_plot = sampling_ids[, 2],
                                                const = list(sampling_points = sampling_points,
                                                             input_layer = nlcd_layer,
                                                             what = "landscape",
                                                             classes_max = 3),
                                                n_jobs = nrow(sampling_ids),
                                                template = list(queue = "medium",
                                                                walltime = "24:00:00",
                                                                mem_cpu = "12G",
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

#### FUTURE (clip and calc) ####

# load the packages
# library("future")
# library("future.batchtools")
# library("furrr")
# 
# # now we specify a future topology that fits our HPC
# # login node -> cluster nodes -> core/ multiple cores
# login <- tweak(future::remote,
#                workers = "gwdu103.gwdg.de",
#                user = "hesselbarth3") # user = login credential
# 
# sbatch <- tweak(future.batchtools::batchtools_slurm,
#                 template = "future_slurm.tmpl",
#                 resources = list(job.name = "calculate_lsm", # name of the job
#                                  log.file = "calculate_lsm.log", # name of log file
#                                  queue = "medium", # which partition
#                                  service = "normal", # which QOS
#                                  walltime = "06:00:00", # walltime <hh:mm:ss>
#                                  processes = 1)) # number of cores
# 
# future::plan(list(login, sbatch, future::sequential)) # how to run on nodes, could also be sequential
# 
# # no max size of globals
# options(future.globals.maxSize = Inf)
# 
# sampling_ids <- suppoRt::expand_grid_unique(x = seq_along(sampling_points), 
#                                             y = seq_along(sampling_points))
#                                            
# landscape_metrics %<-% furrr::future_map2(sampling_ids[, 1], sampling_ids[, 2], 
#                                           function(x, y) {
#                                             clip_and_calc(focal_plot = x, 
#                                                           other_plot = y, 
#                                                           sampling_points = sampling_points,
#                                                           input_layer = nlcd_layer, 
#                                                           what = "landscape", 
#                                                           classes_max = 3)})
# 
# future::resolved(future::futureOf(landscape_metrics))