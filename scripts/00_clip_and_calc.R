clip_and_calc <- function(focal_plot, other_plot, 
                          sampling_points, 
                          input_layer = NULL, path = NULL,
                          ...) {
  
  options(to_disk = TRUE)
  
  # load input_layer on HPC if null
  if (is.null(input_layer)) {
    
    # check if path is specified
    if (is.null(path)) {
      stop("Please provide path to input_layer.rds.", call. = FALSE)
    }
    
    # load layer
    input_layer <- readRDS(file = path) 
  }
  
  # calculate distance matrix
  distance_matrix <- as.matrix(dist(sampling_points@coords, method = "euclidean", 
                                    diag = TRUE, upper = TRUE)) 
  
  dist <- distance_matrix[focal_plot, other_plot]
  
  # center of ellipse between point focal_plot and other_plot
  center_x <- (sampling_points@coords[focal_plot, 1] + 
                 sampling_points@coords[other_plot, 1]) / 2 #Ellipse center x
  
  center_y <- (sampling_points@coords[focal_plot, 2] + 
                 sampling_points@coords[other_plot, 2]) / 2 #Ellipse center y
  
  # angle of ellipsoid
  ellipsoid_angle  <- atan((sampling_points@coords[other_plot, 2] - 
                              sampling_points@coords[focal_plot, 2]) / 
                             (sampling_points@coords[other_plot, 1] - 
                                sampling_points@coords[focal_plot, 1])) # ellipse angle
  
  # distance between the sampling points below 4000 
  # which is the maxium distance the model is trained for
  if (dist < 4000) { 
    
    # use the trained model to predict the radii using the acutal distance between
    # the two samplings points as explanatory variable
    radius_a <- (exp(-0.001163 * distance_matrix[focal_plot,other_plot] + 
                       5.366365) + 60) / 100 * distance_matrix[focal_plot, other_plot]
    
    radius_b <- (exp(-0.001163 * distance_matrix[focal_plot,other_plot] + 
                       5.600736) + 40) / 100 * distance_matrix[focal_plot, other_plot]
  }
  
  # distance larger than model trained for
  if (dist >= 4000) {
    
    radius_a <- 60 / 100 * distance_matrix[focal_plot, other_plot]
    
    radius_b <- 40 / 100 * distance_matrix[focal_plot, other_plot]
  }
  
  # calculate ellipsoid using the model prediction
  radian <- seq(0,2 * pi, length.out = 360)
  
  x_coords <- center_x + radius_a * cos(radian) * 
    cos(ellipsoid_angle) - radius_b * sin(radian) * sin(ellipsoid_angle)
  
  y_coords <- center_y + radius_a * cos(radian) * 
    sin(ellipsoid_angle) + radius_b * sin(radian) * cos(ellipsoid_angle)
  
  # create sp
  ellipsoid <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(x_coords, 
                                                                            y_coords))), 
                                                     ID = 1)))
  
  # clip to ellipsoid
  input_layer  <- raster::crop(input_layer, ellipsoid)
  
  # crop to extent of ellipsoid
  input_layer <- raster::mask(input_layer, ellipsoid) 
  
  result <- landscapemetrics::calculate_lsm(input_layer, progress = FALSE, verbose = FALSE, 
                                            ...)
  
  # add site information
  result <- dplyr::mutate(result, 
                          site_a = focal_plot, 
                          site_b = other_plot, 
                          euclidean_distance = dist)
  
  # result$site_a <- focal_plot
  # 
  # result$site_b <- other_plot
  # 
  # result$euclidean_distance <- dist
  
  return(result)
}
