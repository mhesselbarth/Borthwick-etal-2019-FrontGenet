library(checkmate)
library(raster)
library(SDMTools)
library(tidyverse)

# Function to calculate metrics for each patch within class (used later)
Patch.Labeling <- function(input, cellsize){
  input_1 <- input # Raster for class 1
  input_2 <- input # Raster for class 0
  input_2[] <- abs(input_2[] - 1) # Invert values
  
  input_list <- list(input_1, input_2) # Save as list
  
  input_list %>%
    purrr::map_dfr(., function(input){
      input %>%
        SDMTools::ConnCompLabel() %>% # Assign neighbouring cells of the same class as one patch
        SDMTools::PatchStat(cellsize = cellsize) %>% # Calculate patch metrics
        as.data.frame() %>% # Convert to dataframe
        tibble::as.tibble() %>% # Convert to tibble
        dplyr::filter(patchID!=0)},  # Remove first patch (All cells not belonging to current class)
      .id='class')
}

# Function to convert raster to polygon. The function requires 'gdal_polygonize.py' from GDAL.
# For information: https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
polygonizer <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile', 
                        pypath=NULL, readpoly=TRUE, quietish=TRUE) {
  # x: an R Raster layer, or the file path to a raster file recognised by GDAL
  # outshape: the path to the output shapefile (if NULL, a temporary file will be created)
  # gdalformat: the desired OGR vector format
  # pypath: the path to gdal_polygonize.py (if NULL, an attempt will be made to determine the location
  # readpoly: should the polygon shapefile be read back into R, and returned by this function? (logical)
  # quietish: should (some) messages be suppressed? (logical)
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  ## The line below has been commented:
  # if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.") 
  owd <- getwd()
  on.exit(setwd(owd))
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists)) 
      stop(sprintf('File already exists: %s', 
                   toString(paste(outshape, c('shp', 'shx', 'dbf'), 
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.asc')})
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  
  system2('C:\\OSGeo4W64\\OSGeo4W.bat', 
          args=(sprintf('"%s" "%s" %s -f "ESRI Shapefile" "%s.shp"', 
                        pypath, rastpath, ifelse(quietish, '-q ', ''), outshape)))
  
  
  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quietish)
    return(shp) 
  }
  return(NULL)
}


# Function to calculate the edge length. 
Edge.Length <- function(input){
  input %>%
    polygonizer() %>% 
    sf::st_as_sfc() %>%
    sf::st_cast('MULTILINESTRING') %>% 
    sf::st_union() %>% 
    sf::st_length() %>%
    as.numeric()
}

LandscapeStat <- function(input, cellsize=1){
  
  # Convert input data to list depending on class of input
  if(class(input)=="RasterLayer"){ # input is raster layer
    cat("\nConverted RasterLayer to list\n\n") # write information on console
    input <- raster::as.list(input) # conver to list
  }
  
  else if(class(input)=="RasterStack"){ # input is raster stack (several raster layers)
    cat("\nConverted RasterStack to list\n\n") # write information on console
    input <- raster::as.list(input) # convert to lust
  }

  else{
    if(checkmate::test_list(input)){ # input is alreay list
      cat("\nList provided as input\n\n") # write information on console
    }
    else{stop("\nPlease provide list containing landscape(s) as input\n\n")} # input is none of the above classes - stop function
  }
  
 result_class <- input %>% 
   purrr::map_dfr(function(x){tibble::as.tibble(SDMTools::ClassStat(x, cellsize=cellsize))}, .id="landscape") %>% # calculate class-level metrics
   dplyr::mutate(landscape=as.integer(landscape)) %>% # Clearer name for ID variable 
   dplyr::group_by(landscape) %>% # Grouping for later summarising of results for each landscape
   dplyr::mutate(Proportion_squared = (total.area / sum(total.area)) ^ 2) %>% # Proportion of each class in landscape
   dplyr::summarise(AI = sum(aggregation.index * prop.landscape) / sum(prop.landscape), # Calculate landscape metrics based on classes
                    SIDI = 1 - sum(Proportion_squared),
                    SIEI = dplyr::if_else(condition = is.na(SIDI / (1 - (1 / n()))),
                                          true = 0, 
                                          false = SIDI / (1 - (1 / n()))), 
                    TA = sum(total.area) / 10000) 
   
 result_patch <- input %>%
   purrr::map_dfr(Patch.Labeling, cellsize=cellsize, .id='landscape') %>% # Calculate patch levels metrics
   dplyr::mutate(landscape=as.integer(landscape)) %>% # Clearer name for ID 
   dplyr::group_by(landscape) %>% # Grouping for later summarising of results for each landscape
   dplyr::mutate(perim.area.ratio.w = (perimeter / area) * area) %>%
   dplyr::summarise(AREA_CV = cv(area, aszero=TRUE), # Calculate landscape metrics based on single patches
                    SHAPE_CV = cv(shape.index, aszero=TRUE),
                    FRAC_CV = cv(frac.dim.index, aszero=TRUE),
                    PARA_AM = sum(perim.area.ratio.w) / (sum(area) / 10000),
                    CORE_CV = cv(core.area, aszero=TRUE), 
                    LPI = max(area) / sum(area) * 100, 
                    DIVISION = 1 - sum((area / sum(area)) ^ 2 ),
                    SPLIT = sum(area) / (sum(area ^ 2) / sum(area)))
 
 result_landscape <- input %>%
   purrr::map_dfr(function(x){tibble(TE = Edge.Length(x))}, .id = 'landscape') %>% 
   dplyr::mutate(landscape = as.integer(landscape))
   
 result <- dplyr::full_join(x=result_class,
                                 y=result_patch, 
                                 by='landscape') %>%
   dplyr::full_join(result_landscape, 
                    by = 'landscape') %>%
   dplyr::mutate(ED = TE / TA, 
                 LSI = (0.25 * TE) / sqrt((TA * 10000))) %>%
   dplyr::select(-TE)
 
 return(result) # return result
}
