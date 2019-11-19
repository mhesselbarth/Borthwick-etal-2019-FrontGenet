################################
# Check spatial autocorrelation 
################################

## Install corMLPE from github
devtools::install_github("nspope/corMLPE")

## Load corMLPE package
library(corMLPE)

## Load distances dataframe
df <- read.csv("/Volumes/Hard_Ext/DGS presentations/data/Output/Data_S1.csv", stringsAsFactors=F)
head(df)

## Sort data
df <- rst[order(rst$site_1, rst$site_2),]

## Identify which samples are from same location based on pairwise geographic distance
ulab <- unique(c(df$site_1, df$site_2))
dis <- matrix(0, length(ulab), length(ulab))
rownames(dis) <- colnames(dis) <- ulab
for(i in 1:nrow(df)){
  dis[df$site_1[i],df$site_2[i]] <- dis[df$site_2[i],df$site_1[i]] <- df$euclidean_distance[i]
  location <- cutree(hclust(as.dist(dis)),h=0) #assign individuals to unique locations
}

### Model IBD using MLPE and NMLPE
m1 <- gls(RST ~ euclidean_distance, correlation = corMLPE(form = ~site_1+site_2), data = df) ## MLPE model

acf(resid(m1,type='normalized')) ## No spatial dependence pattern

summary(m1)
