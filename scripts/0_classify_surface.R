# Import libraries
# library(BAMMtools)
library(patchwork) # devtools::install_github("thomasp85/patchwork")
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(raster)
library(sp)
library(tidyverse)

# Import data
habitat_surface <- getwd() %>% 
  paste0("/data/GIS/habitat_surface.tif") %>%
  raster::raster()

# Get all values in raster
raster_values <- na.omit(raster::getValues(habitat_surface))

# Classifiy breaks using Fisher-Jenks algorithm
range(raster_values)
# breaks <- BAMMtools::getJenksBreaks(raster_values, k = 4)

# breaks <- c(-4.687700, -1.798261, 6.580765)                                 # 2 classes
breaks <- c(-4.6876998, -2.7545104, 0.2864172, 6.5807652)                     # 3 classes
# breaks <- c(-4.687700, -3.188930, -1.102762, 0.156716, 1.721571, 6.580765)  # 5 classes

# Cut the values into classes
raster_values_df <- tibble::tibble(value = raster_values)

# For some reason the highest value is NA, I guess because of rounding issues
raster_values_df <- dplyr::mutate(raster_values_df, class = cut(raster_values_df$value,
                                                                breaks = breaks,
                                                                labels = c("low", "medium", "high"),
                                                                include.lowest = TRUE,
                                                                dig.lab = max(nchar(breaks))))
        
# Plot histogramm
ggplot_histogramm <- ggplot2::ggplot(data = raster_values_df, 
                                     aes(x = value)) + 
  ggplot2::geom_histogram(ggplot2::aes(y = ..count.., 
                                       fill = class), 
                          binwidth = 0.25, col = "white") +
  ggplot2::scale_fill_viridis_d(name = "Class") +
  ggplot2::labs(x = "Habitat surface quality", y = "Count") +
  ggplot2::theme_bw()

# Classify original habitats
habitat_surface_pmm <- raster::cut(habitat_surface, 
                                   breaks = breaks, 
                                   include.lowest = TRUE, 
                                   dig.lab = max(nchar(breaks)))

# Plot habitat surface and PMM
ggplot_surface <- ggplot2::ggplot(data = raster::as.data.frame(habitat_surface,
                                                               na.rm = TRUE, 
                                                               xy = TRUE)) + # plot data
    ggplot2::geom_raster(aes(x = x, y = y, fill = habitat_surface)) + # plot patches with fill according to value
    ggplot2::coord_equal() + # ratio of 1:1 of x and y axis
    ggplot2::scale_fill_viridis_c(name = "Habitat quality") + # setting colours of fill
    ggplot2::labs(x = "x coordinate", y = "y-coordinate", title = "Habitat surface") + # set labels axis
    ggplot2::theme_classic() + # overall appearance of the plot
    ggplot2::theme(axis.text=element_blank(), # remove xy-coordinates from axis (not very informative)
                 axis.ticks=element_blank())


ggplot_pmm <- ggplot2::ggplot(data = raster::as.data.frame(habitat_surface_pmm,
                                                           na.rm = TRUE, 
                                                           xy = TRUE)) + # plot data
  ggplot2::geom_raster(aes(x = x, y = y, fill = as.factor(layer))) + # plot patches with fill according to value
  ggplot2::coord_equal() + # ratio of 1:1 of x and y axis
  ggplot2::scale_fill_viridis_d(name = "Habitat quality", labels = c("low", "medium", "high")) + # setting colours of fill
  ggplot2::labs(x = "x coordinate", y = "y-coordinate", title = "PMM") + # set labels axis
  ggplot2::theme_classic() + # overall appearance of the plot
  ggplot2::theme(axis.text=element_blank(), # remove xy-coordinates from axis (not very informative)
                 axis.ticks=element_blank())

ggplot_surface_vs_pmm <- ggplot_surface + ggplot_pmm

# Save output
overwrite <- FALSE

helpeR::save_rds(object = habitat_surface_pmm, 
                 filename = "habitat_surface_pmm.rds", 
                 path = paste0(getwd(), "/data/output"), 
                 overwrite = overwrite)

helpeR::save_ggplot(plot = ggplot_histogramm, 
                    filename = "ggplot_histogramm.png", 
                    path = paste0(getwd(), "/plots"), 
                    overwrite = overwrite)

helpeR::save_ggplot(plot = ggplot_surface_vs_pmm, 
                    filename = "ggplot_surface_vs_pmm.png", 
                    path = paste0(getwd(), "/plots"), 
                    overwrite = overwrite)
