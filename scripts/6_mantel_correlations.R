# load packages
library(ecodist)
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions")
library(tidyverse)


# import data
landscape_metrics <- readr::read_rds(paste0(getwd(), "/data/output/landscape_metrics.rds"))
surface_metrics <- readr::read_rds(paste0(getwd(), "/data/output/surface_metrics.rds"))
rst <- readr::read_rds(paste0(getwd(), "/data/output/rst.rds"))

# only landscape level metrics and NA metrics
landscape_metrics <- dplyr::filter(landscape_metrics, 
                                   level == "landscape", 
                                   !is.na(value))

# add coloumn to identifiy metrics
# landscape_metrics$metric_type <- "patch_metric"
# surface_metrics$metric_type <- "surface_metric"

#reshape surface_metrics to long format
surface_metrics <- tidyr::gather(surface_metrics, 
                                 -name, -landscape, -site_1, -site_2,
                                 key = "metric", value = "value")


# adding RST
landscape_metrics <- left_join(x = landscape_metrics, 
                               y = rst,
                               by = c("site_a" = "site_1", 
                                      "site_b" = "site_2"),
                               suffix = c(".surface_metrics", ".rst"))

surface_metrics <- left_join(x = surface_metrics, 
                             y = rst,
                             by = c("site_1","site_2"), 
                             suffix = c(".surface_metrics", ".rst")) 

# Calculate Mantel correlations
mantel_correlations_surface <- dplyr::group_by(surface_metrics, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = ecodist::mantel(RST ~ value)[[1]])

# adding label for later plotting
mantel_correlations_surface <- dplyr::mutate(mantel_correlations_surface, 
              label = paste0(metric, " - Mantel=", round(mantel_cor, 3)))

mantel_correlations_patch <-dplyr::group_by(landscape_metrics, metric) %>% 
  dplyr::summarise(n = n(),
                   mantel_cor = ecodist::mantel(RST ~ value)[[1]])

surface_metrics <- dplyr::left_join(x = surface_metrics, 
                                    y =mantel_correlations_surface, 
                                    by = "metric")




# mantel correlation
labeller_surface <- function(variable,value){
  return(mantel_correlations_surface$label)
}

ggplot_correlation_rst_surface <- ggplot2::ggplot(surface_metrics, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1) + 
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(~ metric, scales = "free_x", labeller = labeller_surface) + 
  ggplot2::labs(x = "Surface metrics", y = "Rst") + 
  ggplot2::theme_bw()

ggplot_correlation_rst_patch <- ggplot2::ggplot(landscape_metrics, ggplot2::aes(x = value, y = RST)) + 
  ggplot2::geom_point(pch = 1, na.rm = TRUE) + 
  ggplot2::geom_smooth(method = "lm", na.rm = TRUE) +
  ggplot2::facet_wrap(~ metric, scales = "free_x") + 
  ggplot2::labs(x = "Patch metrics", y = "Rst") + 
  ggplot2::theme_bw()
