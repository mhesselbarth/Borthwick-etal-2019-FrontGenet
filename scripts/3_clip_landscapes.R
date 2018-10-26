##############################################################
# Clipping landscape using CRW
# 1. CRW simulation
# 2. Extract ellipse information from success CRWs, then make equation
# 3. Clipping landscape using the ellipse
##############################################################

library(purrr)      # mapping
library(maptools)   # read shape files
library(PBSmapping) # make polygons
library(raster)     # read raster files and clipping
library(spatstat)   # for random walker simulation

##############################################################
#  1. CRW simulation
#   1-1. define 3 functions for CRW simulation
#   1-2. CRW simulation (takes about 2 hours), then draw ellipse (potential success paths)
##############################################################

##############################################################
# 1-1.  define 3 functions for CRW simulation
##############################################################

# Semi-minor radius of ellipse can be obtained by calculating the maximum 
#     perpendicular distance from the line linking between start point 
#     of source patch and hit point of target patch
sminorP2L <- function (x0, y0, x1, y1, xr, yr) {
  
  vx0 <- rep(x0,length(xr))
  vy0 <- rep(y0,length(xr))
  vx1 <- rep(x1,length(xr))
  vy1 <- rep(y1,length(xr))
  
  a = (vy1-vy0)/(vx1-vx0)
  b = vy0 - a*vx0
  if (x1 == x0) distance = abs(xr - vx0)
  if (y1 == y0) distance = abs(yr- vy0)
  if (x1 != x0 & y1 != y0) distance = abs(a*xr -yr + b)/ sqrt(a^2 + 1)
  
  return (max(distance))
}

# Function 2: smajorP2L ()
# Semi-major radius of ellipse can be obtained by calculating the maximum 
#     perpendicular distance from the line seperating between start point 
#     of source patch and hit point of target patch
smajorP2L <- function (x0, y0, x1, y1, xr, yr) {
  
  vx0 <- rep(x0,length(xr))
  vy0 <- rep(y0,length(xr))
  vx1 <- rep(x1,length(xr))
  vy1 <- rep(y1,length(xr))
  
  a = -(vx1-vx0)/(vy1-vy0)
  b = (vy0+vy1)/2 - a*(vx0+vx1)/2
  if (x1 == x0) distance =  abs(yr - (vy0 + vy1)/2)
  if (y1 == y0) distance =  abs(xr - (vx0 + vx1)/2)
  if (x1 != x0 & y1 != y0) distance = abs(a * xr - yr + b)/ sqrt(a^2 + 1)
  
  return (max(distance))
}


# isOnlineseg function find which hitting point is the first one
#  because crossing.psp can't tell this information (just show hitting points)
#  find that random walker's successful path from start to hitting point on target area

isOnlineseg  <- function (x0, y0, x1,y1,hit.px, hit.py) {
  
  hit.steps <- rep(NA, length(hit.px))
  
  for( hit.n in 1:length(hit.px)) {
    
    hit.x   <-  rep(hit.px[hit.n], stepCount)
    hit.y   <-  rep(hit.py[hit.n], stepCount)
    
    hit.steps[hit.n] <- min(
      which(
        (y1 - y0) / (x1 - x0) * (hit.x - x1) + (y1 - hit.y) < 10^(-5)
        & ifelse(x1 > x0, x0 < hit.x & hit.x < x1, x1 < hit.x & hit.x < x0)
        & ifelse(y1 > y0, y0 < hit.y & hit.y < y1, y1 < hit.y & hit.y < y0)
        )
      )
  }
  
  # first line segment hitting target and hitting point vector index in hit.point
  c(min(hit.steps), which(hit.steps==min(hit.steps)))
}

################################################################################
#  1-2. CRW simulation based on different distance between source and target
################################################################################

# The source is always at at point source(0,0)
xs <- 0              # source point x location
ys <- 0              # source point y location

# The target point is at increasing distance in y-direction (but same x-coordinate)
yt <- c(250, 500, 750, 1000, 1500, 2000, 3000, 4000) # target y location
xt <- rep(0, length(yt))                      # target x location

# Plot basic "experiment set-up" - The 'S' is always the starting point of the random walk
# and differente distance to the target 'T' are used (increasing y-coordinate)
ggplot2::ggplot() +
  ggplot2::geom_point(ggplot2::aes(x = xt, y = yt), shape = 'T', size = 2.5) +
  ggplot2::geom_label(ggplot2::aes(x = xt, y = yt, label = yt), hjust = 1.25) +
  ggplot2::geom_point(ggplot2::aes(x = xs, y = ys), shape = 'S', size = 2.5) +
  ggplot2::theme_classic()
  
# Preallocate vector for results for all possible distances (dimension of ellipsoids)
A <-rep(NA, length(yt))        # Eliipse major radius axis
B <-rep(NA, length(yt))        # Eliipse minor radius axis
Eangle <- rep(NA, length(yt))  # Ellipse angle

success.perc <- rep(NA,length(yt)) # store success.percentage information for each target and possible distance

# Define global variavles
stepCount <- 1000    # maximum step count
n.walker  <- 100000  # number of random walker

k <- 0.85            # degree of correlation between movement directions: highly correlated

# plot results next to each other
par(mfrow = c(2,4))

for (target in 1:length(yt)) { # loop for all different distances between source and target
  print(yt[target]) # print current target
  
  set.seed(100000) # set seed for reproducibility
  
  n.success <- 0 # set counter for successful walks for distance to zero 
  
  frame.e <- rep(NA,n.walker) # Preallocate vector for random steps 
  s.ellipse <- data.frame("Ex.c" = frame.e, "Ey.c" = frame.e, "E.angle" = frame.e,
                          "minor.r" = frame.e, "major.r" = frame.e) # Preallocate dataframe for ellipsoids
  
  plot(0, 0, xlim = c(-500, 500), ylim = c(-200, max(yt) + 1000), asp = 1, xlab = "x", ylab = "y") # plot starting point
  text(0, max(yt) + 1000, paste("distance = ",yt[target], " m"), cex = 1.5) # write distance to target
  
  w.extent <- c(-500000, 500000, -500000, 500000)
  
  # strart point: circle polygon
  start.C <-  as.psp(disc(2, c(xs, ys)), w.extent) # The random walk can start at any point on the circle (and not just the point itself)
  start.P  <- runifpointOnLines(n.walker, start.C) # Random points on the starting circle as possible starting points
  
  # target point: circle polygon 
  target.C <- as.psp(disc(2, c(xt[target],yt[target])), w.extent  ) # the source is not just one point but rather a circle the random walk can hit
  
  # add the current source and targe to the plot
  plot(start.C, add = TRUE)
  plot(target.C, add = TRUE)
  
  for (walker in 1:n.walker) { # loop over random walks 
    
    # Create random directions for each step of random walk using a gaussian distribution
    turningAngles <- round(rnorm(stepCount, mean = 0, sd = (1-k) * 2 * pi), 2) # sd means direction of animal head
    turningAngles[1] <- runif(1, 0, 2 * pi) # make sure that first step goes in random direction
    
    # Create random length for each step of random walk using a gamma distribution
    stepLength <- round(rgamma(stepCount, shape = 2, scale = 25), 2) # each step has step length
    
    theta <- cumsum(turningAngles) # theta are now turning angles relative to north
    
    # Multiply the length of each step by the direction of each step in x and y direction
    dx <- stepLength * sin(theta)
    dy <- stepLength * cos(theta)
    
    # Add the random steps to the coordinates of each step before by using the cumulativ sum
    # starting from a point on the target circle
    x <- c(start.P$x[walker], start.P$x[walker] + cumsum(dx))  # now step x,y from start and end
    y <- c(start.P$y[walker], start.P$y[walker] + cumsum(dy))
    
    From <- as.ppp(cbind(x[1:stepCount],y[1:stepCount]), w.extent )
    To <-   as.ppp(cbind(x[2:(stepCount+1)],y[2:(stepCount+1)]), w.extent )
    r.path <- as.psp(from = From, to = To)  # r.path is segement of line xy: psp objet
    
    # Finds if the line of the random walk crosses (hits) the line of the target circle
    (hit.point <- crossing.psp (r.path, target.C)) # to find out hitting point
    
    # when random walker hit the target
    
    if(hit.point$n > 0) { # if the random walk hit the target there will be a crossing, i.e. n>0
      
      # The random walk could hit the targert circle at several points, 
      # but we just interested in the first time it happens
      hit.step1 <- isOnlineseg(r.path$ends$x0, r.path$ends$y0, r.path$ends$x1, r.path$ends$y1,
                               hit.point$x, hit.point$y) # the first hitting point: hit.step1
      
      # Extract the coordinates of the random walk until it hits the source circle for the first time
      success.x <-c(r.path$ends$x0[1:hit.step1[1]], hit.point$x[hit.step1[2]]) # success.x means the r.path from start to hitting step
      success.y <-c(r.path$ends$y0[1:hit.step1[1]], hit.point$y[hit.step1[2]])
      
      # Add line of walk to the graph
      lines(success.x, success.y, col = sample(rainbow(100),1)) # success r.path
      
      # Add the actual start and source point on the circles to the plot
      points(success.x[1],success.y[1])                       # the first step of success r.path
      points(success.x[length(success.x)],success.y[length(success.y)]) # the end step of success r.path
      
      # memorize ellipse information from the success r.path
      # each row represents a random walk, if not successful the row is NA
      s.ellipse[walker,] <- c((success.x[1]+success.x[length(success.x)]) / 2, # Center of ellipsoid in y direction
                               
                               (success.y[1]+success.y[length(success.y)]) / 2, # Center of ellipsoid in y direction
                               
                               atan((success.y[length(success.y)]-success.y[1]) # Angle of ellipsoid
                                    /(success.x[length(success.x)]-success.x[1])),
                               # Calculate minor and major axes of ellipsoid 
                               # See Koh et al. 2013 Fig. 2 for a nice visualization
                               sminorP2L(success.x[1], success.y[1], # Minor axis of ellipsoid
                                         success.x[length(success.x)], success.y[length(success.y)],
                                         success.x,success.y),
                               smajorP2L(success.x[1], success.y[1],
                                         success.x[length(success.x)], success.y[length(success.y)],
                                         success.x,success.y)
      )
      
      # Count sucessful random walk
      n.success <- n.success + 1
      
    } # the end for success hitting random walker
    
  }    # the end for the each random walker simulation
  
  
  # drawing ellipse for success r.paths
  # Using the mean of all successful path to draw ellipsoid containg most of
  # the successful random walks
  Xc <- mean(s.ellipse$Ex.c, na.rm = TRUE)              # mean center point of ellipse
  Yc <- mean(s.ellipse$Ey.c, na.rm = TRUE)              # mean center point of ellipse
  A[target] <- mean(s.ellipse$major.r, na.rm = TRUE)    # A is the major radius
  B[target] <- mean(s.ellipse$minor.r, na.rm = TRUE)    # B is the minor radius
  Eangle[target] <- atan((yt[target]-Yc) / (xt[target] - Xc))# Eliipse angle
  
  t <- seq(0, 2 * pi, 0.1)
  Xe <- Xc + A[target] * cos(t) * cos(Eangle[target]) - B[target] * sin(t) * sin(Eangle[target])
  Ye <- Yc + A[target] * cos(t) * sin(Eangle[target]) + B[target] * sin(t) * cos(Eangle[target])
  
  # Plotting ellipsoid
  lines(Xe, Ye, lty = 2, lwd = 2)
  
  # Calculate rate of successful walks
  success.perc[target] <- n.success / n.walker * 100
  
  # Print results on console
  text(0, max(yt), paste("success =", round(success.perc[target],2), "%"), cex=1.5)
  
} # CRW simulation end.


###############################################################################
# 2. Extract ellipse information from success CRWs, then make equation
#  2-1. ploting the relationship between distance and ellipse raidus.
#  2-2. make equation from the relation
###############################################################################

###############################################################################
#  2-1. plotting the relationship between distance and ellipse raidus
###############################################################################
#distance <- dist.st[which(!is.na(A))]

# All possible distances
distance <- yt

# Sucess-rates of different distances
success.p <- success.perc[1:length(distance)]

# Plot rate against distance
plot(distance, success.p, ylab="Success rate of 100,000 random walker (%)")

# Relative radius of ellpisoid in relation to distance between points 
major.r <- A/distance * 100
minor.r <- B/distance * 100

# Plot results
# See Koh et al. 2013 Fig. 3 for another example
plot(distance, major.r,  ylim = c(0,max(major.r, minor.r)),
     pch = 16,cex = 1.5,
     ylab = "Length of radius in proportion to the distance (%)",
     xlab = "Distance between source and target")
points(distance, minor.r, pch = 17,cex = 1.5)
abline(h = 50, lty = 2)
legend("topright",
       legend = c("major radius (A)", "minor radius (B)"),
       pch = c(16,17))

###############################################################################
#  2-2. make equation from the relation
##############################################################################

# We assume minor and major radius converge to 40% and 60

# Model explainig the minor radius of all sucessful random walks with the distance as explanatory variabel
model.minor <- lm(log(minor.r - 40) ~ distance)

# Using the model to predict radii for all distance between 1:10000
y.minor <- exp(predict(model.minor,list(distance = 1:10000))) + 40
lines(1:10000, y.minor, lty = 2) # add line to plot

# same for major radius
model.major <- lm(log(major.r - 60) ~ distance )
y.minor <- exp(predict(model.major,list(distance = 1:10000))) + 60
lines(1:10000, y.minor, lty = 3)

# store model equation for later us
# minor radius  =  (exp(-0.001163*distance + 5.600736) + 40) * distance (between sourse and target)
# major radius  =  (exp(-0.001163*distance + 5.366365) + 60) * distance



################################################################################
# 3. Clipping landscape with the ellipse
##################################################################################

# Input data for clipping landscape. This is where we have to throw in 
# our patch data forest vs non-forest

# input_layer <- raster(paste0(getwd(), "/data/GIS/forest_300m_0_1.tif"))
# input_layer[values(input_layer) < -999] <- NA
input_layer <- readRDS(paste0(getwd(), "/data/output/habitat_surface_pmm.rds"))

plot(input_layer)

# Read sampling points
#spoints <- readShapeSpatial(paste0(path, "/GIS_data/Sampling_points/sampling_points.shp"))
spoints <- raster::shapefile(paste0(getwd(), "/data/GIS/SSR_17_sites.shp"))

## important !! id: 0 -> 70 
id <- as.numeric(rownames(spoints@data)) # setting unique id of sampling points

plot(spoints, add = TRUE, pch = 1) # plot

# creat distance matrix 
# between all sample points -> diagonal is 0 because thats the same point
# everything above and below the diagonal is equal cause 1 <-> 2 == 2 <-> 1
dist.m <- as.matrix(dist(spoints@coords, method = "euclidean", 
                         diag = TRUE, upper = TRUE)) 

# Loop for all sampling points
clippings <- purrr::map(1:length(id), function(focal_plot) {

  cat(paste0("\rProgress: ", focal_plot , " from ", length(id)))
    
  non_focal_ids <- which(id != focal_plot & id > focal_plot)
  
  # loop for next point (other_plot) to current point (focal_plot) 
  result <- purrr::map(non_focal_ids, function(other_plot) {
    
    # Center of ellipse between point focal_plot and other_plot
    Xc <- (spoints@coords[focal_plot,1] + spoints@coords[other_plot,1]) / 2 #Ellipse center x
    Yc <- (spoints@coords[focal_plot,2] + spoints@coords[other_plot,2]) / 2 #Ellipse center y
    
    # Angle of ellipsoid
    Eangle  <- atan(
      (spoints@coords[other_plot,2] - spoints@coords[focal_plot,2]) / 
        (spoints@coords[other_plot,1] - spoints@coords[focal_plot,1])
    ) # ellipse angle
    
    # distance between the sampling points below 4000 
    # which is the maxium distance the model is trained for
    if(dist.m[focal_plot,other_plot] < 4000) { 
      
      # Use the trained model to predict the radii using the acutal distance between
      # the two samplings points as explanatory variable
      A <- (exp(-0.001163 * dist.m[focal_plot,other_plot] + 5.366365) + 60)/100 * 
        dist.m[focal_plot,other_plot]
      B <- (exp(-0.001163 * dist.m[focal_plot,other_plot] + 5.600736) + 40)/100 * 
        dist.m[focal_plot,other_plot]
    }
    
    # Distance larger than model trained for
    if(dist.m[focal_plot,other_plot] >= 4000) {
      A <- 60/100 * dist.m[focal_plot,other_plot]
      B <- 40/100 * dist.m[focal_plot,other_plot]
    }
    
    # Calculate ellipsoid using the model prediction
    radian <- seq(0,2 * pi, length.out = 360)
    Xe <- Xc + A * cos(radian) * cos(Eangle) - B * sin(radian)*sin(Eangle)
    Ye <- Yc + A * cos(radian) * sin(Eangle) + B * sin(radian)*cos(Eangle)
    Eline <- cbind(Xe, Ye)
    
    # make ellipse polygon from line by two steps
    Epolyset <- as.PolySet(data.frame(PID = rep(1,length(radian)),
                                      SID = rep(1,length(radian)), POS = 1:length(radian),
                                      X = Eline[,"Xe"], Y = Eline[,"Ye"]),
                           projection = 1)
    Epolygon <- PolySet2SpatialPolygons(Epolyset)
    
    # Plot
    # plot(Epolygon, add = TRUE)
    
    # Clipping Ellipse
    # Emask   <- rasterize(Epolygon, input_layer)
    Elands  <- crop(input_layer, Epolygon)
    Elands <- mask(Elands, Epolygon) # Crop to extend of Emask
    
    names(Elands) <- paste0("comparison_", focal_plot, "_", other_plot)
    
    return(Elands)
    
    # Write files
    # writeRaster(Elands, filename=paste("Elands_",focal_plot,"_",other_plot,sep=""),format="GTiff", overwrite=TRUE)
  })
  
  return(result)
})

clippings_flatten <- purrr::flatten(clippings)

UtilityFunctions::save_rds(clippings_flatten, filename = "clippings_pmm.rds", 
                           path = paste0(getwd(), "/data/output"), 
                           overwrite = TRUE)

