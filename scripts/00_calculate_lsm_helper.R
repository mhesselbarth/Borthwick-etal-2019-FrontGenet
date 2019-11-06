# Wrapper around calculate_lsm that not all options must be set for Q()

calculate_lsm_helper <- function(landscape, what, classes_max) {
  
  landscapemetrics::calculate_lsm(landscape = landscape, 
                                  what = what,
                                  classes_max = classes_max,
                                  verbose = FALSE, 
                                  progress = FALSE)
}
