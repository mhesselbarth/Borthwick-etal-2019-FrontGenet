## ----message=FALSE, warning=TRUE-----------------------------------------
require(LandGenCourse)
require(GeNetIt)
require(raster)
require(gdistance)
require(rgdal)

## ----Read in suitability data--------------------------------------------------------------------
#data(rasters, package="GeNetIt")
#RasterMaps <- raster::stack(rasters)
hsl<-raster("c:/users/richard/downloads/T.v_habitat_surface.tif")
Raster<-stack(hsl)

## ------read in sample location data------------------------------------------------------------------
sites<-readOGR("c:/users/richard/downloads/SSR_17_sites.shp")

#reading in Hossam's conductance values:
cond<-raster("c:/users/richard/downloads/t.v_conduct.asc")
#for interest
plot(cond)
plot(hsl)
par(mfrow=c(1,2))

## ----fig.width=8, fig.height=5.5-----------------------------------------
raster::plot(Raster)

## ------------------------------------------------------------------------
raster::plot(Raster$T.v_habitat_surface)
points(sites, pch=3)

## ------------------------------------------------------------------------
#cti <- raster::resample(cti, gsp, method= "bilinear")

## ------------------------------------------------------------------------
#RasterMaps$err27

## ------------------------------------------------------------------------
#err.cost <- (1/RasterMaps$err27)
#err.cost
##Need to revalue everything to get NA as 0 and neg values as positive.

Raster$T.v_habitat_surface<-Raster$T.v_habitat_surface+5
cost<-1/Raster$T.v_habitat_surface
cst.cond<-1/cond$T.v_conduct

## ------------------------------------------------------------------------
#RasterMaps$ffp

## ------------------------------------------------------------------------
#ffp.cost <- (RasterMaps$ffp/5)
#ffp.cost

## ------------------------------------------------------------------------
#RasterMaps$gsp

## ------------------------------------------------------------------------
#gsp.cost <- (RasterMaps$gsp-196)/15
#gsp.cost

## ------------------------------------------------------------------------
#RasterMaps$cti

## ------------------------------------------------------------------------
#cti.cost <- RasterMaps$cti/5
#cti.cost

## ------------------------------------------------------------------------
#cost1 <- (gsp.cost + cti.cost + err.cost + ffp.cost)
#cost1

## ------------------------------------------------------------------------
tr.cost1 <- gdistance::transition(cost, transitionFunction=mean, directions=8) 
tr.cost.cond<-gdistance::transition(cst.cond, transitionFunction=mean, directions=8)

## ------------------------------------------------------------------------
par(mar=c(2,2,1,1))
raster::plot(raster::raster(tr.cost1))
raster::plot(raster::raster(tr.cost.cond))

## ------------------------------------------------------------------------
tr.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)
tr.cost.cond <- gdistance::geoCorrection(tr.cost.cond,type = "c",multpl=FALSE)

## ------------------------------------------------------------------------
par(mar=c(2,2,1,2))
AtoB <- gdistance::shortestPath(tr.cost1, origin=sites[1,], 
                                goal=sites[2,], output="SpatialLines")
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)",legend.lab="Conductance")
lines(AtoB, col="red", lwd=2)
points(sites[1:2,])

par(mar=c(2,2,1,2))
AtoB <- gdistance::shortestPath(tr.cost.cond, origin=sites[1,], 
                                goal=sites[2,], output="SpatialLines")
raster::plot(raster::raster(tr.cost.cond), xlab="x coordinate (m)", 
             ylab="y coordinate (m)",legend.lab="Conductance")
lines(AtoB, col="red", lwd=2)
points(sites[1:2,])

## ----message=FALSE-------------------------------------------------------
par(mar=c(2,2,1,2))
raster::plot(raster::raster(tr.cost1), xlab="x coordinate (m)", 
             ylab="y coordinate (m)", legend.lab="Conductance")
points(sites)

Neighbours <- spdep::tri2nb(sites@coords, row.names = sites$SiteName)

plot(Neighbours, sites@coords, col="darkgrey", add=TRUE)
for(i in 1:length(Neighbours))
{
  for(j in Neighbours[[i]][Neighbours[[i]] > i])
  {
    AtoB <- gdistance::shortestPath(tr.cost1, origin=sites[i,], 
                                goal=sites[j,], output="SpatialLines")
    lines(AtoB, col="red", lwd=1.5)
  }
}

par(mar=c(2,2,1,2))
raster::plot(raster::raster(tr.cost.cond), xlab="x coordinate (m)", 
             ylab="y coordinate (m)", legend.lab="Conductance")
points(sites)

Neighbours <- spdep::tri2nb(sites@coords, row.names = sites$SiteName)

plot(Neighbours, sites@coords, col="darkgrey", add=TRUE)
for(i in 1:length(Neighbours))
{
  for(j in Neighbours[[i]][Neighbours[[i]] > i])
  {
    AtoB <- gdistance::shortestPath(tr.cost.cond, origin=sites[i,], 
                                    goal=sites[j,], output="SpatialLines")
    lines(AtoB, col="red", lwd=1.5)
  }
}

## ------------------------------------------------------------------------
cost1.dist <- gdistance::costDistance(tr.cost1,sites)
cond.dist<-gdistance::costDistance(tr.cost.cond,sites)

## ------------------------------------------------------------------------
comm1.dist <- gdistance::commuteDistance(x = tr.cost1, coords = sites)
comm.cond.dist<-gdistance::commuteDistance(tr.cost.cond,sites)
## ------------------------------------------------------------------------
dist_df <- data.frame("cost1.dist"=as.numeric(cost1.dist),
                      "comm1.dist"=as.numeric(comm1.dist))
dist_cond<- data.frame("cost1.dist"=as.numeric(cond.dist),
                       "comm1.dist"=as.numeric(comm.cond.dist))
## ------------------------------------------------------------------------
corr.LCD.comm <- cor(dist_df$cost1.dist, dist_df$comm1.dist, method = "spearman")
corr.LCD.comm
plot(cost1.dist~comm1.dist)

corr.LCD.cond <- cor(dist_cond$cost1.dist, dist_cond$comm1.dist, method = "spearman")
corr.LCD.cond
plot(cond.dist~comm.cond.dist)

## ------------------------------------------------------------------------
cor_cost <- c()
cor_comm <- c()
res_fact <- seq(2,20,2)
res_fact2<- seq(2,22,2)
for(fac in res_fact){
  cost1_agg <- raster::aggregate(cost, fact = fac)
  tr.cost_agg <- gdistance::transition(cost1_agg, 
                 transitionFunction=mean, directions=8)
  tr.cost_agg <- gdistance::geoCorrection(tr.cost_agg,type = "c",multpl=FALSE)
  cost.dist_agg <- gdistance::costDistance(tr.cost_agg,sites)
  comm.dist_agg <- gdistance::commuteDistance(x = tr.cost_agg, coords = sites)
  cost.dist_agg <- as.numeric(cost.dist_agg)
  comm.dist_agg <- as.numeric(comm.dist_agg)
  cor_cost <- c(cor_cost,cor(dist_df$cost1.dist, cost.dist_agg, 
                             method = "spearman"))
  cor_comm <- c(cor_comm,cor(dist_df$comm1.dist, comm.dist_agg, 
                             method = "spearman"))
}
#conductance
for(fac in res_fact2){
  cost2_agg <- raster::aggregate(cst.cond, fact = fac)
  tr.cost2_agg <- gdistance::transition(cost2_agg, 
                                       transitionFunction=mean, directions=8)
  tr.cost2_agg <- gdistance::geoCorrection(tr.cost2_agg,type = "c",multpl=FALSE)
  cost.dist2_agg <- gdistance::costDistance(tr.cost2_agg,sites)
  comm.dist2_agg <- gdistance::commuteDistance(x = tr.cost2_agg, coords = sites)
  cost.dist2_agg <- as.numeric(cost.dist2_agg)
  comm.dist2_agg <- as.numeric(comm.dist2_agg)
  cor_cost2 <- c(cor_cost,cor(dist_cond$cost1.dist, cost.dist2_agg, 
                             method = "spearman"))
  cor_comm2 <- c(cor_comm,cor(dist_cond$comm1.dist, comm.dist2_agg, 
                             method = "spearman"))
}
## ------------------------------------------------------------------------
par(mar=c(4,4,1,1))
plot(y = cor_cost, x = res_fact, col = "red", pch = 19, 
     ylim = c(0,1), xlab = "Aggregation factor", ylab = "Spearman correlation")
points(y = cor_comm, x = res_fact, col = "blue", pch = 19)
legend("bottomleft", legend = c("Costdist","Commdist"), 
       pch = 19, col = c("red", "blue"))

par(mar=c(4,4,1,1))
plot(y = cor_cost2, x = res_fact2, col = "red", pch = 19, 
     ylim = c(0,1), xlab = "Aggregation factor", ylab = "Spearman correlation")
points(y = cor_comm2, x = res_fact2, col = "blue", pch = 19)
legend("bottomleft", legend = c("Costdist","Commdist"), 
       pch = 19, col = c("red", "blue"))


## ----message=FALSE, warning=TRUE, include=FALSE--------------------------
#dist_df is the data, now need to link it to RST:
rst <- read.csv("C:/Users/richard/Downloads/TV_gene_flow.csv")
View(rst)
rst<-rst[,1:2]
#make a flow/cost df
fc_df<-cbind(dist_df,rst)
cond_df<-cbind(dist_cond,rst)
#check correlations:
lmod.001<-lm(fc_df$RST~fc_df$cost1.dist)
summary(lmod.001)
plot(fc_df$cost1.dist,fc_df$RST)

lmod.002<-lm(fc_df$RST~fc_df$comm1.dist)
summary(lmod.002)
plot(fc_df$comm1.dist,fc_df$RST)

lmod.003<-lm(cond_df$RST~cond_df$cost1.dist)
summary(lmod.003)
plot(cond_df$cost1.dist,cond_df$RST)

lmod.004<-lm(cond_df$RST~cond_df$comm1.dist)
summary(lmod.004)
plot(cond_df$comm1.dist,cond_df$RST)

##doing same with shard alleles
paf <- read.csv("C:/Users/richard/Downloads/Prop.Sh.Allele.csv")
View(paf)
cond_df<-cbind(cond_df,paf$Prop.SH.Al)

lmod.005<-lm(cond_df$`paf$Prop.SH.Al`~cond_df$cost1.dist)
summary(lmod.005)
plot(cond_df$`paf$Prop.SH.Al`,cond_df$cost1.dist)

lmod.006<-lm(cond_df$`paf$Prop.SH.Al`~cond_df$comm1.dist)
summary(lmod.006)
plot(cond_df$comm1.dist,cond_df$`paf$Prop.SH.Al`)

##make matrices:
flow<-dist(fc_df$RST)
cst<-dist(fc_df$cost1.dist)
comm<-dist(fc_df$comm1.dist)


flow2<-as.matrix(flow)
comm2<-as.matrix(comm)
cst2<-as.matrix(cst)

flow.c<-dist(cond_df$RST)
cst.c<-dist(cond_df$cost1.dist)
comm.c<-dist(cond_df$comm1.dist)


flow3<-as.matrix(flow.c)
comm3<-as.matrix(comm.c)
cst3<-as.matrix(cst.c)
##Mantel:
library(ncf)
partial.mantel.test(flow2,comm2,cst2,resamp=1000,method="pearson")
partial.mantel.test(flow3,comm3,cst3,resamp=1000,method="pearson")

