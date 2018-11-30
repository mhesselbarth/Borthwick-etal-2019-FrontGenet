# Wrapper around calculate_lsm that not all options must be set for Q()

calculate_lsm_helper <- function(landscape, level, classes_max) {
  landscapemetrics::calculate_lsm(landscape = landscape, 
                                  level = level,
                                  classes_max = classes_max,
                                  verbose = FALSE)
}
