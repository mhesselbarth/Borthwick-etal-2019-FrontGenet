# load packages
library(ecodist)
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions")
library(tidyverse)

# import data
landscape_metrics <- readr::read_rds(paste0(getwd(), "/data/output/landscape_metrics.rds"))
surface_metrics <- readr::read_rds(paste0(getwd(), "/data/output/surface_metrics.rds"))
rst <- readr::read_rds(paste0(getwd(), "/data/output/rst.rds"))

# only landscape level metrics that are comparable between landscape difference size and no NA values
# absolute_metrics <- c("ca", "ndca", "np", "pafrac", "pr", "ta", "tca", "te")

landscape_level <- dplyr::filter(landscape_metrics,
                                 level == "landscape",
                                 !is.na(value))

# class_level_forest <- dplyr::filter(landscape_metrics,
#                                     level == "class",
#                                     class == 1,
#                                     !is.na(value))
# 
# class_level_complementary <- dplyr::filter(landscape_metrics,
#                                            level == "class",
#                                            class == 2,
#                                            !is.na(value))


#reshape surface_metrics to long format
surface_metrics <- tidyr::gather(surface_metrics, 
                                 -name, -landscape, -site_1, -site_2,
                                 key = "metric", value = "value")


# adding RST
landscape_level <- left_join(x = landscape_level, 
                             y = rst,
                             by = c("site_a" = "site_1", 
                                    "site_b" = "site_2"),
                             suffix = c(".surface_metrics", ".rst"))

surface_metrics <- left_join(x = surface_metrics, 
                             y = rst,
                             by = c("site_1","site_2"), 
                             suffix = c(".surface_metrics", ".rst")) 

# return NA in case ecodist::mantel() fails
safe_mantel <- purrr::possibly(ecodist::mantel, otherwise = NA)

# calculate Mantel correlations
mantel_correlations_surface <- dplyr::group_by(surface_metrics, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = safe_mantel(RST ~ value)[[1]])

# adding label for later plotting
mantel_correlations_surface <- dplyr::mutate(mantel_correlations_surface, 
              label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

# calculate Mantel correlations
mantel_correlations_patch <- dplyr::group_by(landscape_level, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = safe_mantel(RST ~ value)[[1]])

# adding label for later plotting
mantel_correlations_patch <- dplyr::mutate(mantel_correlations_patch, 
                                           label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

# helper function to plot results
labeller_surface <- function(variable,value){
  return(mantel_correlations_surface$label)
}

labeller_patch <- function(variable,value){
  return(mantel_correlations_patch$label)
}

# plot results
ggplot_correlation_rst_surface <- ggplot2::ggplot(surface_metrics, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1) + 
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_surface) + 
  ggplot2::labs(x = "Surface metrics", y = "Rst") + 
  ggplot2::theme_bw()

ggplot_correlation_rst_patch <- ggplot2::ggplot(landscape_level, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1, na.rm = TRUE) + 
  ggplot2::geom_smooth(method = "lm", na.rm = TRUE) +
  ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_patch) + 
  ggplot2::labs(x = "Patch metrics", y = "Rst") + 
  ggplot2::theme_bw()

# save plots
overwrite <- FALSE

UtilityFunctions::save_ggplot(plot = ggplot_correlation_rst_surface, 
                              filename = "ggplot_correlation_rst_surface.png", 
                              path = paste0(getwd(), "/plots"),
                              overwrite = overwrite)

UtilityFunctions::save_ggplot(plot = ggplot_correlation_rst_patch, 
                              filename = "ggplot_correlation_rst_patch.png", 
                              path = paste0(getwd(), "/plots"),
                              overwrite = overwrite)
