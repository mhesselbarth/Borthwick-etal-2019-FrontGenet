# Wrapper around calculate_lsm that not all options must be set for Q()

calculate_lsm_helper <- function(landscape, what) {
  landscapemetrics::calculate_lsm(landscape = landscape, 
                                  what = what,
                                  verbose = FALSE)
}
