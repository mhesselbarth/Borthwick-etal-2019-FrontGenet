# load libraries
library(clustermq)
library(landscapemetrics)
library(raster)
library(tidyverse)

source(paste0(getwd(), "/scripts/0_calculate_lsm_helper.R"))

# load the clippings
clippings_pmm <- readRDS(paste0(getwd(), "/data/output/clippings_pmm_nlcd.rds"))

# check if all rasters all loaded in memory
all(purrr::map_lgl(clippings_pmm, raster::inMemory))

# extract names
names_clippings <- purrr::map_chr(clippings_pmm, function(x) names(x))
names_clippings <- stringr::str_split(names_clippings, pattern = "_", simplify = TRUE) # need for local version

class <- c("lsm_c_ai", 
           "lsm_c_area_mn", 
           "lsm_c_ca", 
           "lsm_c_cai_mn", 
           "lsm_c_clumpy", 
           "lsm_c_core_mn", 
           "lsm_c_cpland", 
           "lsm_c_division", 
           "lsm_c_iji", 
           "lsm_c_lpi", 
           "lsm_c_lsi", 
           "lsm_c_mesh", 
           "lsm_c_np", 
           "lsm_c_pd", 
           "lsm_c_pladj", 
           "lsm_c_pland", 
           "lsm_c_split", 
           "lsm_c_te")

landscape <- c("lsm_l_ai", 
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
               "lsm_l_split", 
               "lsm_l_ta", 
               "lsm_l_te")

what <- c(class, landscape)

# # Calculate metrics locally
# landscape_metrics <- landscapemetrics::calculate_lsm(clippings_pmm,
#                                                      what = what,
#                                                      classes_max = 3)
# 
# # Calculate metrics locally but overall printing progress
# landscape_metrics <- purrr::map(seq_along(clippings_pmm), function(x) {
# 
#   print(paste0("Progress: ", x, " from ", length(clippings_pmm)))
# 
#   calculate_lsm(landscape = clippings_pmm[[x]],
#                 what = what,
#                 classes_max = 3,
#                 verbose = FALSE,
#                 progress = FALSE)
# })
# 
# # Add name of sites
# landscape_metrics <- dplyr::mutate(landscape_metrics,
#                                    site_a = as.integer(names_clippings[layer, 2]),
#                                    site_b = as.integer(names_clippings[layer, 3]))

# Calculate metrics on high performance cluster
landscape_metrics <- clustermq::Q(fun = calculate_lsm_helper,
                                  landscape = clippings_pmm,
                                  const = list(what = what,
                                               classes_max = 3),
                                  n_jobs = length(clippings_pmm),
                                  template = list(queue = "mpi",
                                                  walltime = "48:00",
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

# # FUTURE for HPC
# # load the packages
# library("future")
# library("future.batchtools")
# library("furrr")
# 
# # now we specify a future topology that fits our HPC
# # login node -> cluster nodes -> core/ multiple cores
# login <- future::tweak(remote, workers = "gwdu103.gwdg.de",
#                        user = "hesselbarth3")
# 
# bsub <- future::tweak(batchtools_lsf, template = "lsf.tmpl",
#                       resources = list(job.name = "lsm_clippings",
#                                        log.file = "lsm_clippings.log",
#                                        queue = "mpi",
#                                        walltime = "48:00",
#                                        processes = 1)) # hpc -> nodes
# 
# future::plan(list(login, bsub, future::sequential)) # how to run on nodes, could also be sequential
# 
# landscape_metrics %<-% furrr::future_map(clippings_pmm, function(x) {
#   
#   calculate_lsm(landscape = x,
#                 what = what,
#                 classes_max = 3,
#                 verbose = FALSE,
#                 progress = FALSE)
# })
# 
# future::resolved(future::futureOf(landscape_metrics))