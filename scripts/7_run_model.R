# load libraries
library(tidyverse)
library(lme4)

# import data surface metrics
surface_metrics <- readr::read_rds("data/Output/surface_metrics.rds")

# import data landscape metrics
landscape_metrics <- readr::read_rds("data/Output/landscape_metrics.rds")

# only metrics on landscape level and needed cols
landscape_metrics_lndscp <- dplyr::filter(landscape_metrics, 
                                          level == "landscape") %>%
  dplyr::select(site_a, site_b, metric, value, euclidean_distance)

# reshape to wide format
landscape_metrics_lndscp <- tidyr::spread(landscape_metrics_lndscp, 
                                          metric, value)
