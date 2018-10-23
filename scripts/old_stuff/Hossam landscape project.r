##############################################################
# Clipping landscape using CRW
# 1. CRW simulation
# 2. Extract ellipse information from success CRWs, then make equation
# 3. Clipping landscape using the ellipse
##############################################################


library(spatstat)   # for random walker simulation
library(maptools)   # read shape files
library(raster)     # read raster files and clipping
library(PBSmapping) # make polygons


##############################################################
#  1. CRW simulation
#   1-1. define 3 functions for CRW simulation
#   1-2. CRW simulation (takes about 2 hours), then draw ellipse (potential success paths)
##############################################################

##############################################################
# 1-1.  define 3 functions for CRW simulation
##############################################################

# calculate max minor axis distance for success random walker
# x0,y0 are source point, x1,y1 are target point
sminorP2L <- function (x0, y0, x1, y1, xr, yr) {

          vx0 <- rep(x0,length(xr))
          vy0 <- rep(y0,length(xr))
          vx1 <- rep(x1,length(xr))
          vy1 <- rep(y1,length(xr))

          a = (vy1-vy0)/(vx1-vx0)
          b = vy0 - a*vx0
          if (x1==x0) distance = abs(xr - vx0)
          if (y1==y0) distance = abs(yr- vy0)
          if (x1!=x0 & y1!=y0) distance = abs(a*xr -yr + b)/ sqrt(a^2 + 1)

          return (max(distance))
          }

# calculate max major axis distance for success random walker
smajorP2L <- function (x0, y0, x1, y1, xr, yr) {

          vx0 <- rep(x0,length(xr))
          vy0 <- rep(y0,length(xr))
          vx1 <- rep(x1,length(xr))
          vy1 <- rep(y1,length(xr))

          a = -(vx1-vx0)/(vy1-vy0)
          b = (vy0+vy1)/2 - a*(vx0+vx1)/2
          if (x1==x0) distance =  abs(yr - (vy0+vy1)/2)
          if (y1==y0) distance =  abs(xr - (vx0+vx1)/2)
          if (x1!=x0 & y1!=y0) distance = abs(a*xr -yr + b)/ sqrt(a^2 + 1)

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

                  hit.steps[hit.n] <-
                  min (which((y1- y0)/(x1-x0)*(hit.x - x1) + (y1 - hit.y) < 10^(-5)
                                     & ifelse(x1>x0, x0 < hit.x & hit.x < x1, x1 < hit.x & hit.x < x0)
                                     & ifelse(y1>y0, y0 < hit.y & hit.y < y1, y1 < hit.y & hit.y < y0)
                                     )
                      )
                  }

                  # first line segment hitting target and hitting point vector index in hit.point
                  c(min(hit.steps), which(hit.steps==min(hit.steps)))

                  }

################################################################################
#  1-2. CRW simulation based on different distance between source and target
################################################################################

xs <- 0              # source point x location
ys <- 0              # source point y location

xt <- rep(0, length(yt))                      # target x location
yt <- c(250,500,750,1000,1500,2000,3000,4000) # target y location


A <-rep(NA, length(yt))        # Eliipse major radius axis
B <-rep(NA, length(yt))        # Eliipse minor radius axis
Eangle <- rep(NA, length(yt))  # Ellipse angle

success.perc <- rep(NA,length(yt)) # store success.percentage information for each target


stepCount <- 1000    # maximum step count
n.walker  <- 100000  # number of random walker

k <- 0.85            # degree of correlation between movement directions: highly correlated


par(mfrow=c(2,4))

for ( target in 1: length(yt)) {
 print(yt[target])

 set.seed(100000)

 n.success <- 0

 frame.e <- rep(NA,n.walker)
 s.ellipse <- data.frame("Ex.c"=frame.e, "Ey.c"=frame.e, "E.angle"=frame.e,
                          "minor.r"=frame.e, "major.r"=frame.e) # success.ellipse

 plot(0,0,xlim=c(-500,500), ylim=c(-200,max(yt)+1000), asp=1,xlab="x", ylab="y")
 text(0, max(yt)+1000, paste("distance = ",yt[target], " m"), cex=1.5)

  w.extent <- c(-500000,500000,-500000,500000)

  # strart point: circle polygon
  start.C <-  as.psp(disc(2, c(xs,ys)), w.extent)
  start.P  <- runifpointOnLines(n.walker, start.C)

  # target point: circle polygon
  target.C <- as.psp(disc(2, c(xt[target],yt[target])), w.extent  )

  plot(start.C, add=T)
  plot(target.C, add=T)

  for (walker in 1:n.walker) {

    turningAngles <- round(rnorm(stepCount, mean=0, sd=(1-k)*2*pi),2) # sd means direction of animal head
    turningAngles[1] <- runif(1, 0, 2*pi) # make sure that first step goes in random direction

    stepLength <- round(rgamma(stepCount, shape=2, scale=25),2) # each step has step length

    theta <- cumsum(turningAngles) # theta are now turning angles relative to north

    dx <- stepLength * sin(theta)
    dy <- stepLength * cos(theta)

    x <- c(start.P$x[walker], start.P$x[walker] + cumsum(dx))  # now step x,y from start and end
    y <- c(start.P$y[walker], start.P$y[walker] + cumsum(dy))

    From <- as.ppp(cbind(x[1:stepCount],y[1:stepCount]), w.extent )
    To <-   as.ppp(cbind(x[2:(stepCount+1)],y[2:(stepCount+1)]), w.extent )
    r.path <- as.psp(from=From, to=To)  # r.path is segement of line xy: psp objet

    hit.point <- crossing.psp (r.path, target.C) # to find out hitting point

    # when random walker hit the target
    
    if(hit.point$n > 0) {
    
        hit.step1 <- isOnlineseg(r.path$ends$x0, r.path$ends$y0, r.path$ends$x1, r.path$ends$y1,
                   hit.point$x, hit.point$y) # the first hitting point: hit.step1

        success.x <-c(r.path$ends$x0[1:hit.step1[1]],hit.point$x[hit.step1[2]]) # success.x means the r.path from start to hitting step
        success.y <-c(r.path$ends$y0[1:hit.step1[1]],hit.point$y[hit.step1[2]])

        lines(success.x, success.y, col=sample(rainbow(100),1)) # success r.path
        points(success.x[1],success.y[1])                       # the first step of success r.path
        points(success.x[length(success.x)],success.y[length(success.y)]) # the end step of success r.path

        # memorize ellipse information from the success r.path
        s.ellipse[walker,] <- c( (success.x[1]+success.x[length(success.x)])/2,

                                (success.y[1]+success.y[length(success.y)])/2,

                                atan((success.y[length(success.y)]-success.y[1])
                                      /(success.x[length(success.x)]-success.x[1])),

                                sminorP2L (success.x[1], success.y[1],
                                           success.x[length(success.x)], success.y[length(success.y)],
                                           success.x,success.y),
                                smajorP2L (success.x[1], success.y[1],
                                           success.x[length(success.x)], success.y[length(success.y)],
                                           success.x,success.y)
                               )

        n.success <- n.success + 1

     } # the end for success hitting random walker

  }    # the end for the each random walker simulation


#drawing ellipse for success r.paths
Xc <- mean(s.ellipse$Ex.c, na.rm=T)              # mean center point of ellipse
Yc <- mean(s.ellipse$Ey.c, na.rm=T)              # mean center point of ellipse
A[target] <- mean(s.ellipse$major.r, na.rm=T)    # A is the major radius
B[target] <- mean(s.ellipse$minor.r, na.rm=T)    # B is the major radius
Eangle[target] <- atan((yt[target]-Yc)/(xt[target]-Xc))# Eliipse angle

t <- seq(0,2*pi,0.1)
Xe <- Xc + A[target]*cos(t)*cos(Eangle[target]) - B[target]*sin(t)*sin(Eangle[target])
Ye <- Yc + A[target]*cos(t)*sin(Eangle[target]) + B[target]*sin(t)*cos(Eangle[target])

lines(Xe, Ye, lty=2, lwd=2)


success.perc[target] <- n.success / n.walker * 100

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
distance <- yt

success.p <- success.perc[1:length(distance)]
plot(distance, success.p, ylab="Success rate of 100,000 random walker (%)")

major.r <- A/distance * 100
minor.r <- B/distance * 100

plot(distance, major.r,  ylim=c(0,max(major.r, minor.r)),
      pch=16,cex=1.5,
      ylab="Length of radius in proportion to the distance (%)",
      xlab="Distance between source and target")
points(distance, minor.r, pch=17,cex=1.5)
abline(h=50, lty=2)
legend("topright",
       legend=c("major radius (A)", "minor radius (B)"),
        pch=c(16,17))
        
###############################################################################
#  2-2. make equation from the relation
##############################################################################

# We assume minor and major radius converge to 40% and 60
model.minor <- lm(log(minor.r - 40) ~ distance )
y.minor<-exp(predict(model.minor,list(distance=1:10000)))+40
lines(1:10000, y.minor, lty=2)


model.major <- lm(log(major.r - 60) ~ distance )
y.minor<-exp(predict(model.major,list(distance=1:10000)))+60
lines(1:10000, y.minor, lty=3)

# store model equation

# minor radius  =  (exp(-0.001163*distance + 5.600736) + 40) * distance (between sourse and target)
# major radius  =  (exp(-0.001163*distance + 5.366365) + 60) * distance



################################################################################
# 3. Clipping landscape with the ellipse
##################################################################################

setwd("C:/Hossam_surface/")

Indi <- raster("T_surface.tif")
plot(Indi)

spoints <- readShapeSpatial(paste0(getwd(), "/Sampling_points/sampling_points.shp"))

## important !! id: 0 -> 70
id <- as.numeric(rownames(spoints@data))

plot(spoints, add=T, pch=1)

# creat distance matrix
dist.m <- as.matrix(dist(spoints@coords, method="euclidean", diag=T, upper=T))

for (i in 1:length(id)) {

  for (j in (i+1):length(id)) {

        Xc <- (spoints@coords[i,1] + spoints@coords[j,1])/2 #Ellipse center x
        Yc <- (spoints@coords[i,2] + spoints@coords[j,2])/2 #Ellipse center x

        Eangle  <- atan(
                  (spoints@coords[j,2]-spoints@coords[i,2])/(spoints@coords[j,1]-spoints@coords[i,1])
                   ) # ellipse angle

        if(dist.m[i,j] < 4000) {

            A <- (exp(-0.001163*dist.m[i,j] + 5.366365) + 60)/100 * dist.m[i,j]
            B <- (exp(-0.001163*dist.m[i,j] + 5.600736) + 40)/100 * dist.m[i,j]
        }
        if(dist.m[i,j] >= 4000) {
        
            A <- 60/100 * dist.m[i,j]
            B <- 40/100 * dist.m[i,j]
        }
        
        radian <- seq(0,2*pi, length.out=360)
        Xe <- Xc + A*cos(radian)*cos(Eangle) - B*sin(radian)*sin(Eangle)
        Ye <- Yc + A*cos(radian)*sin(Eangle) + B*sin(radian)*cos(Eangle)
        Eline <- cbind(Xe, Ye)

        # make ellipse polygon from line by two steps
        Epolyset <- as.PolySet(data.frame(PID=rep(1,length(radian)),
                        SID=rep(1,length(radian)), POS=1:length(radian),
                        X= Eline[,"Xe"], Y= Eline[,"Ye"]),
                        projection=1)
        Epolygon <- PolySet2SpatialPolygons(Epolyset)
        plot(Epolygon, add=T)
         
        # Clipping Ellipse
        Emask   <- rasterize(Epolygon, Indi)
        Elands  <- mask(Indi, Emask)
        # Write files
        writeRaster(Elands, filename=paste("Elands_",i,"_",j,sep=""),format="GTiff", overwrite=TRUE)
        
  }

}

