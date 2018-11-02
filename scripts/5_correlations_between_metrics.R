# load packages
library(ggplot2)
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions")
library(tidyverse)

# import data
landscapemetrics <- readr::read_rds(paste0(getwd(), "/data/output/landscape_metrics.rds"))

# Split class level and landscape level
class_level <- dplyr::filter(landscapemetrics, level == "class")
landscape_level <- dplyr::filter(landscapemetrics, level == "landscape")

# plot correlations
# ggplot_correlation_class <- landscapemetrics::show_correlation(class_level)

ggplot_correlation_landscape <- landscapemetrics::show_correlation(landscape_level)

# convert to wide format
metrics_landscape_wide <- stats::xtabs(value ~ layer + metric,
                                       data = landscapemetrics[, c(1, 5:6)])

# attr(metrics_landscape_wide, "landscape") <- NULL
# attr(metrics_landscape_wide, "call") <- NULL

# correlation between metrics
correlation_matrix_landscape <- cor(metrics_landscape_wide, 
                                    method = "pearson")


# set diagonal to  NA
correlation_matrix_landscape[upper.tri(correlation_matrix_landscape,
                                       diag = TRUE)] <- NA

# convert to tibble
correlation_matrix_landscape_df <- tibble::tibble(
  metric_1 = rownames(correlation_matrix_landscape)[row(correlation_matrix_landscape)],
  metric_2 = colnames(correlation_matrix_landscape)[col(correlation_matrix_landscape)],
  value = c(correlation_matrix_landscape)
  )

# remove all NA cases
correlation_landscape <- correlation_matrix_landscape_df[complete.cases(correlation_matrix_landscape_df),]

dplyr::filter(correlation_landscape, value < 0.5)

# Save plots
overwrite <- TRUE

UtilityFunctions::save_ggplot(plot = ggplot_correlation_landscape, 
                              filename = "ggplot_correlation_landscape.png", 
                              path = paste0(getwd(), "/plots"), 
                              width = 35, height = 35, units = "cm",
                              overwrite = overwrite)
