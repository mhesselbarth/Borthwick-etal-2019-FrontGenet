# load packages
library(ecodist)
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(tidyverse)

# import data
landscape_metrics <- readr::read_rds(paste0(getwd(), "/data/output/landscape_metrics.rds"))
surface_metrics <- readr::read_rds(paste0(getwd(), "/data/output/surface_metrics.rds"))
rst <- readr::read_rds(paste0(getwd(), "/data/rst.rds"))

# filter data 
# remove metrics with too high VIF (see 7_run_model.R)
landscape_level <- dplyr::filter(landscape_metrics,
                                 level == "landscape",
                                 !is.na(value), 
                                 metric %in% c("pd", "iji", "prd", "core_mn", 
                                               "cai_mn", "split", "mesh"))

# class_forest <- dplyr::filter(landscape_metrics,
#                               level == "class",
#                               class == 1,
#                               !is.na(value))

# class_complementary <- dplyr::filter(landscape_metrics,
#                                      level == "class",
#                                      class == 2,
#                                      !is.na(value))

# reshape surface_metrics to long format
surface_metrics <- tidyr::gather(surface_metrics, 
                                 -name, -landscape, -site_1, -site_2,
                                 key = "metric", value = "value")

# Remove Sku because of VIF (see 7_run_model.R)
surface_metrics <- dplyr::filter(surface_metrics, 
                                 metric != "Sku")

# adding RST
landscape_level <- left_join(x = landscape_level, 
                             y = rst,
                             by = c("site_a" = "site_1", 
                                    "site_b" = "site_2"))

# class_forest <- left_join(x = class_forest, 
#                           y = rst,
#                           by = c("site_a" = "site_1", 
#                                  "site_b" = "site_2"))

# class_complementary <- left_join(x = class_complementary, 
#                                  y = rst,
#                                  by = c("site_a" = "site_1"))

surface_metrics <- left_join(x = surface_metrics, 
                             y = rst,
                             by = c("site_1","site_2"), 
                             suffix = c(".surface_metrics", ".rst")) 

# return NA in case ecodist::mantel() fails
safe_mantel <- purrr::possibly(ecodist::mantel, otherwise = NA)

# calculate Mantel correlations surface metrics
mantel_correlations_surface <- dplyr::group_by(surface_metrics, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = safe_mantel(RST ~ value)[[1]])

# adding label for later plotting
mantel_correlations_surface <- dplyr::mutate(mantel_correlations_surface, 
              label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

# calculate Mantel correlations landscape level
mantel_correlations_landscape_level <- dplyr::group_by(landscape_level, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = safe_mantel(RST ~ value)[[1]])

# adding label for later plotting
mantel_correlations_landscape_level <- dplyr::mutate(mantel_correlations_landscape_level, 
                                                     label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

# calculate Mantel correlations forest class
# mantel_correlations_forest <- dplyr::group_by(class_forest, metric) %>% 
#   dplyr::summarise(n = n(),
#                    mantel_cor = safe_mantel(RST ~ value)[[1]])

# adding label for later plotting
# mantel_correlations_forest <- dplyr::mutate(mantel_correlations_forest, 
#                                             label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

# # calculate Mantel correlations complementary class
# mantel_correlations_complementary <- dplyr::group_by(class_complementary, metric) %>% 
#   dplyr::summarise(n = n(),
#                    mantel_cor = safe_mantel(RST ~ value)[[1]])
# 
# # adding label for later plotting
# mantel_correlations_complementary <- dplyr::mutate(mantel_correlations_complementary, 
#                                                    label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))



# helper function to plot results
labeller_surface <- function(variable, value){
  return(mantel_correlations_surface$label)
}

labeller_landscape_level <- function(variable, value){
  return(mantel_correlations_landscape_level$label)
}

# labeller_class_forest <- function(variable, value){
#   return(mantel_correlations_forest$label)
# }

# labeller_class_complementary <- function(variable, value){
#   return(mantel_correlations_complementary$label)
# }

# plot results
ggplot_correlation_rst_surface <- ggplot2::ggplot(surface_metrics, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1) + 
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_surface) + 
  ggplot2::labs(x = "Surface metrics", y = "Rst") + 
  ggplot2::theme_bw(base_size = 20)

ggplot_correlation_rst_landscape_level <- ggplot2::ggplot(landscape_level, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1, na.rm = TRUE) + 
  ggplot2::geom_smooth(method = "lm", na.rm = TRUE) +
  ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_landscape_level) + 
  ggplot2::labs(x = "Patch metrics", y = "Rst") + 
  ggplot2::theme_bw(base_size = 20) 

# ggplot_correlation_rst_class_forest <- ggplot2::ggplot(class_forest, ggplot2::aes(x = value, y = RST)) + 
#   ggplot2::geom_point(pch = 1, na.rm = TRUE) + 
#   ggplot2::geom_smooth(method = "lm", na.rm = TRUE) +
#   ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_class_forest) + 
#   ggplot2::labs(x = "Patch metrics", y = "Rst", title = "Class: Forest") + 
#   ggplot2::theme_bw()

# ggplot_correlation_rst_class_complementary <- ggplot2::ggplot(class_forest, ggplot2::aes(x = value, y = RST)) + 
#   ggplot2::geom_point(pch = 1, na.rm = TRUE) + 
#   ggplot2::geom_smooth(method = "lm", na.rm = TRUE) +
#   ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_class_complementary) + 
#   ggplot2::labs(x = "Patch metrics", y = "Rst", title = "Class: Complementary") + 
#   ggplot2::theme_bw()

# save plots
overwrite <- TRUE

helpeR::save_ggplot(plot = ggplot_correlation_rst_surface, 
                    filename = "ggplot_correlation_rst_surface.png", 
                    path = paste0(getwd(), "/plots"),
                    overwrite = overwrite, 
                    width = 15, height = 7, unit = "in")

helpeR::save_ggplot(plot = ggplot_correlation_rst_landscape_level, 
                    filename = "ggplot_correlation_rst_landscape_level.png", 
                    path = paste0(getwd(), "/plots"),
                    overwrite = overwrite,
                    width = 15, height = 7, unit = "in")

# helpeR::save_ggplot(plot = ggplot_correlation_rst_class_forest,
#                     filename = "ggplot_correlation_rst_class_forest.png",
#                     path = paste0(getwd(), "/plots"),
#                     overwrite = overwrite,
#                     width = 15, height = 7, unit = "in")

# helpeR::save_ggplot(plot = ggplot_correlation_rst_class_complementary, 
#                               filename = "ggplot_correlation_rst_class_complementary.png", 
#                               path = paste0(getwd(), "/plots"),
#                               overwrite = overwrite)

