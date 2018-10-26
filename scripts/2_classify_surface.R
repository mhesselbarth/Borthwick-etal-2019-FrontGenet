# Import libraries
library(BAMMtools)
library(patchwork) # devtools::install_github("thomasp85/patchwork")
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions")
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

breaks <- BAMMtools::getJenksBreaks(raster_values, k = 6)

# breaks <- c(-4.687700, -1.798261, 6.580765)                                 # 2 classes
breaks <- c(-4.6876998, -2.7545104, 0.2864172, 6.5807652)                     # 3 classes
# breaks <- c(-4.687700, -3.188930, -1.102762, 0.156716, 1.721571, 6.580765)  # 5 classes

# Plot histogramm
ggplot_histogramm <- ggplot2::ggplot(data = data.frame(value = raster_values), 
                                     aes(x = value)) + 
  ggplot2::geom_histogram(ggplot2::aes(y = ..density..), binwidth = 0.25,
                          fill = "black", col = "white") +
  ggplot2::geom_density(ggplot2::aes(x = value), alpha = 0.65, 
                        fill = "gray", col = "black", size = 1) + 
  ggplot2::geom_vline(xintercept = breaks[1], col = "red") +
  ggplot2::geom_vline(xintercept = breaks[2], col = "red") +
  ggplot2::geom_vline(xintercept = breaks[3], col = "red") +
  ggplot2::geom_vline(xintercept = breaks[4], col = "red") +
  ggplot2::geom_text(aes(x = (breaks[1] + breaks[2]) / 2, y = 0.5, label = "Class 1")) +
  ggplot2::geom_text(aes(x = (breaks[2] + breaks[3]) / 2, y = 0.5, label = "Class 2")) +
  ggplot2::geom_text(aes(x = (breaks[3] + breaks[4]) / 2, y = 0.5, label = "Class 3")) +
  theme_bw()

# Classify original habitats
habitat_surface_pmm <- raster::cut(habitat_surface, 
                                   breaks = breaks, include.lowest = TRUE)

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
  ggplot2::scale_fill_viridis_d(name = "Habitat quality") + # setting colours of fill
  ggplot2::labs(x = "x coordinate", y = "y-coordinate", title = "PMM") + # set labels axis
  ggplot2::theme_classic() + # overall appearance of the plot
  ggplot2::theme(axis.text=element_blank(), # remove xy-coordinates from axis (not very informative)
                 axis.ticks=element_blank())

ggplot_comparison <- ggplot_surface + ggplot_pmm

# Save output
overwrite <- FALSE

saveRDS(object = habitat_surface_pmm, 
        file = paste0(getwd(), "/data/output/habitat_surface_pmm.rds"))

ggplot2::ggsave(plot = ggplot_histogramm, 
                filename = "ggplot_histogramm.png", 
                path = paste0(getwd(), "/data/output"))

ggplot2::ggsave(plot = ggplot_comparison, 
                filename = "ggplot_comparison.png", 
                path = paste0(getwd(), "/data/plots"))
