# load libraries
library(landscapemetrics)
library(raster)
library(sf)
library(tidyverse)

# Full file names of all clipped ellipsoids
clippings_names <- getwd() %>%
  paste0('/data/Clippings') %>%
  list.files(pattern='.tif$', full.names=TRUE)

# Import all tif-files in folder clipping
clippings_list <- clippings_names %>%
  purrr::map(function(x) raster::raster(x))

# Only part of file name containing the pairwise sampling sites
pairwise_sampling_sites <- clippings_names %>%
  purrr::map_chr(function(x){x %>%
      stringr::str_sub(start = stringr::str_locate(., 'Elands_')[2]+1,
                       end = stringr::str_length(.)-4)})

# Only part of file name containing the pairwise sampling sites in form Elands_x_x
pairwise_sampling_sites_Elands <- clippings_names %>%
  purrr::map_chr(function(x){x %>%
      stringr::str_sub(start = stringr::str_locate(., 'Elands_')[1],
                       end = stringr::str_length(.)-4)})

# Calculate landscape-level metrics
landscape_metrics <- clippings_list %>%
  LandscapeStat(cellsize = 300) %>% # Calculate landscape metrics
  dplyr::mutate(landscape = pairwise_sampling_sites, 
                landscape_Elands = pairwise_sampling_sites_Elands) %>% # Add name of sites
  tidyr::separate(landscape, into=c("site_1", "site_2"), sep='_', remove=F) %>% # Split into Site_1 and Site_2
  dplyr::mutate(site_1 = as.integer(site_1), # Convert to integer
                site_2 = as.integer(site_2))


system.time(metrics <- calculate_metrics(clippings_list[[1]], what = "landscape"))

# Import sample points
sample_points <- path %>% 
  paste0("/GIS_data/SSR_17_Sites.shp") %>%
  raster::shapefile() %>%
  as.data.frame() %>%
  tibble::as.tibble()

# Calculate distance between each point
distance_matrix <- sample_points %>%
  dplyr::select(E, N) %>% # Only select the coordinates
  dist(diag = TRUE, upper = TRUE) %>% # Calculate distance matrix
  as.matrix() %>%
  as.data.frame() %>%
  tibble::as.tibble()

# Add euclidean distance to each pair of sitres
landscape_metrics <- landscape_metrics %>%
  dplyr::mutate(euclid_dist=purrr::pmap_dbl(., function(site_1, site_2, ...){ # Apply function to each row of tibble and add coloumn
    distance_matrix[site_1,site_2] %>% # Take site_1 and site_2 of each row and get according value from distance matrix
      as.double()
    }))

anyNA(landscape_metrics) # Check for NAs

landscape_metrics <- landscape_metrics %>%
  dplyr::arrange(site_1, site_2)

readr::write_delim(landscape_metrics, # Write to comma separeted file
                   path=paste0(path, '/landscape_metrics.csv'), 
                   delim = ';')
                 