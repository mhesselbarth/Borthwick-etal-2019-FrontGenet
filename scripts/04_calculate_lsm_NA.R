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
# clippings_pmm_nlcd_NA <- readRDS(paste0(getwd(), "/data/output/clippings_pmm_nlcd_NA.rds"))
# 
# # check if all rasters all loaded in memory
# all(purrr::map_lgl(clippings_pmm_nlcd_NA, raster::inMemory))
# 
# # extract names
# names_clippings <- purrr::map_chr(clippings_pmm_nlcd_NA, function(x) names(x))
# 
# names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

# load input layer
nlcd_layer_NA <- readRDS(paste0(getwd(), "/data/Output/nlcd_reclassified_NA.rds"))

# load sampling points
sampling_points <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))

# #### Calculate locally ####
# 
# # Calculate metrics locally but overall printing progress
# total_clippigings <- length(clippings_pmm_nlcd_NA)
# 
# landscape_metrics <- purrr::map(seq_along(clippings_pmm_nlcd_NA), function(x) {
# 
#   print(paste0("Progress: ", x, " from ", total_clippigings))
# 
#   result <- calculate_lsm(landscape = clippings_pmm_nlcd_NA[[x]],
#                           level = "landscape",
#                           classes_max = 2,
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
landscape_metrics_NA <- suppoRt::submit_to_cluster(fun = clip_and_calc,
                                                   focal_plot = sampling_ids[, 1],
                                                   other_plot = sampling_ids[, 2],
                                                   n_jobs = nrow(sampling_ids),
                                                   log_worker	= TRUE,
                                                   const = list(sampling_points = sampling_points,
                                                                # input_layer = nlcd_layer_NA,
                                                                what = "landscape",
                                                                classes_max = 2, 
                                                                path = "/home/uni08/hesselbarth3/nlcd_reclassified_NA.rds"),
                                                   template = list(job_name = "lsm_clip", 
                                                                   log_file = "lsm_clip.log",
                                                                   queue = "medium",
                                                                   walltime = "48:00:00",
                                                                   mem_cpu = "16384",
                                                                   processes = 1))

suppoRt::save_rds(object = landscape_metrics_NA,
                  filename = "landscape_metrics_NA_raw.rds",
                  path = paste0(getwd(), "/data/Output"),
                  overwrite = FALSE)

# bind to one dataframe 
landscape_metrics_NA <- dplyr::bind_rows(landscape_metrics_NA)

# replace layer with 1:136
landscape_metrics_NA$layer <- rep(x = 1:nrow(sampling_ids), 
                                  each = length(unique(landscape_metrics_NA$metric)))

suppoRt::save_rds(object = landscape_metrics_NA,
                  filename = "landscape_metrics_NA.rds",
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
#                                                           input_layer = nlcd_layer_NA, 
#                                                           what = "landscape", 
#                                                           classes_max = 3)})
# 
# future::resolved(future::futureOf(landscape_metrics))