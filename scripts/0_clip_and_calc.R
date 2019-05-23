clip_and_calc <- function(focal_plot, sampling_points, raster, ...) {
  
  id <- seq_along(sampling_points)
  
  # calculate distance matrix
  distance_matrix <- as.matrix(dist(sampling_points@coords, method = "euclidean", 
                                    diag = TRUE, upper = TRUE)) 
  
  non_focal_ids <- which(id != focal_plot & id > focal_plot)
  
  # loop for next point (other_plot) to current point (focal_plot) 
  result_inner <- purrr::map_dfr(non_focal_ids, function(other_plot) {
    
    dist <- distance_matrix[focal_plot, other_plot]
    
    # Center of ellipse between point focal_plot and other_plot
    center_x <- (sampling_points@coords[focal_plot, 1] + 
                   sampling_points@coords[other_plot, 1]) / 2 #Ellipse center x
    
    center_y <- (sampling_points@coords[focal_plot, 2] + 
                   sampling_points@coords[other_plot, 2]) / 2 #Ellipse center y
    
    # Angle of ellipsoid
    ellipsoid_angle  <- atan((sampling_points@coords[other_plot, 2] - 
                                sampling_points@coords[focal_plot, 2]) / 
                               (sampling_points@coords[other_plot, 1] - 
                                  sampling_points@coords[focal_plot, 1])) # ellipse angle
    
    # distance between the sampling points below 4000 
    # which is the maxium distance the model is trained for
    if (dist < 4000) { 
      
      # Use the trained model to predict the radii using the acutal distance between
      # the two samplings points as explanatory variable
      radius_a <- (exp(-0.001163 * distance_matrix[focal_plot,other_plot] + 
                         5.366365) + 60) / 100 * distance_matrix[focal_plot, other_plot]
      
      radius_b <- (exp(-0.001163 * distance_matrix[focal_plot,other_plot] + 
                         5.600736) + 40) / 100 * distance_matrix[focal_plot, other_plot]
    }
    
    # Distance larger than model trained for
    if (dist >= 4000) {
      
      radius_a <- 60 / 100 * distance_matrix[focal_plot, other_plot]
      
      radius_b <- 40 / 100 * distance_matrix[focal_plot, other_plot]
    }
    
    # Calculate ellipsoid using the model prediction
    radian <- seq(0,2 * pi, length.out = 360)
    
    x_coords <- center_x + radius_a * cos(radian) * 
      cos(ellipsoid_angle) - radius_b * sin(radian) * sin(ellipsoid_angle)
    
    y_coords <- center_y + radius_a * cos(radian) * 
      sin(ellipsoid_angle) + radius_b * sin(radian) * cos(ellipsoid_angle)
    
    # create sp
    ellipsoid <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(x_coords, 
                                                                              y_coords))), 
                                                       ID = 1)))
  
    # Clipping Ellipse
    # Emask   <- rasterize(Epolygon, input_layer)
    ellipsoid_cropped  <- raster::crop(input_layer, ellipsoid)
    
    ellipsoid_cropped <- raster::mask(ellipsoid_cropped, ellipsoid) # Crop to extend of Emask
    
    result_metrics <- landscapemetrics::calculate_lsm(ellipsoid_cropped, ...)
    
    result_metrics <- dplyr::mutate(result_metrics, 
                                    site_a = focal_plot,
                                    site_b = other_plot, 
                                    euclidean_distance = dist)
    
    return(result_metrics)
  })
  
  return(result_inner)
}
