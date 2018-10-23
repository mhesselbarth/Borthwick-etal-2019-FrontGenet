library(ecodist)
require(ecodist)
library(readr)
require(readr)
library(LandGenCourse)
library(ggplot2)
combined_dataset <- read_csv("c:/users/richard/Downloads/combined_datasetV3.csv") 
#there's probably a more efficient way (i.e. make list and iterate over it) but I just added our variables and assigned
#vectors of pairwise distances because I'm rushed for time and not great with R:

#if you want to see a list of the columns in dataset > spec(combined_dataset)

RST<-combined_dataset$RST #genetic distance 

#-----------------------------------------------------------------------------------------
#gradientmetrics and mantel correlation with RST:

#For each of the metrics:
#first create vector
#then run mantel correlation
#then plot variables (RST as y)
#then add linear regression
#then add r value above graph

par(mfrow=c(4,3), oma=c(2,2,1,0)) # mar=c(2,2,2.5,2))


Sa<-combined_dataset$Sa
RST_Sa<-mantel(RST~Sa)
plot(RST~Sa, xlab=NA, ylab=NA)
abline(lm(RST~Sa))
rval<-round(RST_Sa[[1]],3)
mtext(rval,3)
mtext("Sa",1,2.2)
mtext("RST",2,2.2)

S10z<-combined_dataset$S10z
RST_S10z<-mantel(RST~S10z)
plot(RST~S10z, xlab=NA, ylab=NA)
abline(lm(RST~S10z))
rval<-round(RST_S10z[[1]],3)
mtext(rval,3)
mtext("S10z",1,2.2)
mtext("RST",2,2.2)

Ssk<-combined_dataset$Ssk
mantel(RST~Ssk)
RST_Ssk<-plot(RST~Ssk, xlab=NA, ylab=NA)
abline(lm(RST~Ssk))
rval<-round(RST_Ssk[[1]],3)
mtext(rval,3)
mtext("Ssk",1,2.2)
mtext("RST",2,2.2)

Sku<-combined_dataset$Sku
RST_Sku<-mantel(RST~Sku)
plot(RST~Sku, xlab=NA, ylab=NA)
abline(lm(RST~Sku))
rval<-round(RST_Sku[[1]],3)
mtext(rval,3)
mtext("Sku",1,2.2)
mtext("RST",2,2.2)

Sdr<-combined_dataset$Sdr
RST_Sdr<-mantel(RST~Sdr)
plot(RST~Sdr, xlab=NA, ylab=NA)
abline(lm(RST~Sdr))
rval<-round(RST_Sdr[[1]],3)
mtext(rval,3)
mtext("Sdr",1,2.2)
mtext("RST",2,2.2)

Sbi<-combined_dataset$Sbi
RST_Sbi<-mantel(RST~Sbi)
plot(RST~Sbi, xlab=NA, ylab=NA)
abline(lm(RST~Sbi))
rval<-round(RST_Sbi[[1]],3)
mtext(rval,3)
mtext("Sbi",1,2.2)
mtext("RST",2,2.2)

Std<-combined_dataset$Std
RST_Std<-mantel(RST~Std)
plot(RST~Std, xlab=NA, ylab=NA)
abline(lm(RST~Std))
rval<-round(RST_Std[[1]],3)
mtext(rval,3)
mtext("Std",1,2.2)
mtext("RST",2,2.2)

Stdi<-combined_dataset$Stdi
RST_Stdi<-mantel(RST~Stdi)
plot(RST~Stdi, xlab=NA, ylab=NA)
abline(lm(RST~Stdi))
rval<-round(RST_Stdi[[1]],3)
mtext(rval,3)
mtext("Stdi",1,2.2)
mtext("RST",2,2.2)

Sfd<-combined_dataset$Sfd
RST_Sfd<-mantel(RST~Sfd)
plot(RST~Sfd, xlab=NA, ylab=NA)
abline(lm(RST~Sfd))
rval<-round(RST_Sfd[[1]],3)
mtext(rval,3)
mtext("Sfd",1,2.2)
mtext("RST",2,2.2)

Srwi<-combined_dataset$Srwi
RST_Srwi<-mantel(RST~Srwi)
plot(RST~Srwi, xlab=NA, ylab=NA)
abline(lm(RST~Srwi))
rval<-round(RST_Srwi[[1]],3)
mtext(rval,3)
mtext("Srwi",1,2.2)
mtext("RST",2,2.2)
title("Mantel correlation between surface metrics and RST", outer=TRUE)

#-------------------------------------------------------------------------------------
#patchmetrics and mantel correlation with RST:

#For each of the metrics:
#first create vector
#then run mantel correlation
#then plot variables (RST as y)
#then add linear regression
#then add r value above graph

#HAVING SOME ISSUES WITH SOME PATCH METRIC VECTOR MATRICES AS NOT BEING 
#"SQUARE"DUE TO NA VALUES - HOW DO WE FIX THIS? THIS RESULTS IN MANTEL R 
#NOT BEING CALULATED.

par(mfrow=c(4,3), oma=c(2,2,1,1))


AI<-combined_dataset$AI
RST_AI<-mantel(RST~AI)
plot(RST~AI, xlab=NA, ylab=NA)
abline(lm(RST~AI))
rval<-round(RST_AI[[1]],3)
mtext(rval,3)
mtext("AI",1,2.2)
mtext("RST",2,2.2)

SIDI<-combined_dataset$SIDI
RST_SIDI<-mantel(RST~SIDI)
plot(RST~SIDI, xlab=NA, ylab=NA)
abline(lm(RST~SIDI))
rval<-round(RST_SIDI[[1]],3)
mtext(rval,3)
mtext("SIDI",1,2.2)
mtext("RST",2,2.2)

SIEI<-combined_dataset$SIEI
RST_SIEI<-mantel(RST~SIEI)
plot(RST~SIEI, xlab=NA, ylab=NA)
abline(lm(RST~SIEI))
rval<-round(RST_SIEI[[1]],3)
mtext(rval,3)
mtext("SIEI",1,2.2)
mtext("RST",2,2.2)


SHAPE_CV<-combined_dataset$SHAPE_CV
RST_SHAPE_CV<-mantel(RST~SHAPE_CV)
plot(RST~SHAPE_CV, xlab=NA, ylab=NA)
abline(lm(RST~SHAPE_CV))
rval<-round(RST_SHAPE_CV[[1]],3)
mtext(rval,3)
mtext("SHAPE_CV",1,2.2)
mtext("RST",2,2.2)

FRAC_CV<-combined_dataset$FRAC_CV
RST_FRAC_CV<-mantel(RST~FRAC_CV)
plot(RST~FRAC_CV, xlab=NA, ylab=NA)
abline(lm(RST~FRAC_CV))
rval<-round(RST_FRAC_CV[[1]],3)
mtext(rval,3)
mtext("FRAC_CV",1,2.2)
mtext("RST",2,2.2)

PARA_AM<-combined_dataset$PARA_AM
RST_PARA_AM<-mantel(RST~PARA_AM)
plot(RST~PARA_AM, xlab=NA, ylab=NA)
abline(lm(RST~PARA_AM))
rval<-round(RST_PARA_AM[[1]],3)
mtext(rval,3)
mtext("PARA_AM",1,2.2)
mtext("RST",2,2.2)

CORE_CV<-combined_dataset$CORE_CV
RST_CORE_CV<-mantel(RST~CORE_CV)
plot(RST~CORE_CV, xlab=NA, ylab=NA)
abline(lm(RST~CORE_CV))
rval<-round(RST_CORE_CV[[1]],3)
mtext(rval,3)
mtext("CORE_CV",1,2.2)
mtext("RST",2,2.2)

LPI<-combined_dataset$LPI
RST_LPI<-mantel(RST~LPI)
plot(RST~LPI, xlab=NA, ylab=NA)
abline(lm(RST~LPI))
rval<-round(RST_LPI[[1]],3)
mtext(rval,3)
mtext("LPI",1,2.2)
mtext("RST",2,2.2)

DIVISION<-combined_dataset$DIVISION
RST_DIVISION<-mantel(RST~DIVISION)
plot(RST~DIVISION, xlab=NA, ylab=NA)
abline(lm(RST~DIVISION))
rval<-round(RST_DIVISION[[1]],3)
mtext(rval,3)
mtext("DIVISION",1,2.2)
mtext("RST",2,2.2)

SPLIT<-combined_dataset$SPLIT
RST_SPLIT<-mantel(RST~SPLIT)
plot(RST~SPLIT, xlab=NA, ylab=NA)
abline(lm(RST~SPLIT))
rval<-round(RST_SPLIT[[1]],3)
mtext(rval,3)
mtext("SPLIT",1,2.2)
mtext("RST",2,2.2)

title("Mantel correlation between patch metric and RST", outer=TRUE)

#--------------------------------------------------------------------------------------

#Geodistances:
HOSSAM_Geo_dist<-combined_dataset$Geo_dist #Hossamsdata
MAX_euclid_dist<-combined_dataset$euclid_dist #Max's code
ALIDA_DISTANCE<-combined_dataset$ALIDA_DISTANCE #My ArcMap point distance calculation


#Plot geodistances against each other to show that Hossam and our data differ:

par(mfrow=c(1,2), oma=c(2,0,2,0))  #plot two graphs on one page, add margin for title

#plot(Geo_dist~euclid_dist) with simple regression line:
lmgeo1<- lm(HOSSAM_Geo_dist~MAX_euclid_dist) #regression
plot(HOSSAM_Geo_dist~MAX_euclid_dist, xlab=NA, ylab=NA) #actual plot
abline(lmgeo1) #fit regression
rforlmgeo1m <- mantel(HOSSAM_Geo_dist~MAX_euclid_dist)
mtext(rforlmgeo1m[[1]],3) #show r value above plot
mtext("Hossam_vs_Max/Alida_distance",1,2)

plot(MAX_euclid_dist~ALIDA_DISTANCE, xlab=NA, ylab=NA) #didn't even do a regression cause values are identical
mtext("Max_vs_Alida_distance",1,2)
title("Different geographic distances", outer=TRUE)
mtext("rval=1",3)


#--------------------------------------------------------------------------------------

#CALCULATING SPEARMAN RANK CORRELATION BETWEEN PATCH METRICS: ##Richard thought: can we remove all of this and focus on Pearson?

##Richard follow-up thought, we can get a matrix like this: cor.table<-cor(combined_dataset[,19:31],method="spearman") or use 
#method="pearson" then we could ditch all the print text, I think.


# alternative way to calculate Spearman rank correlation (produces same results as mantel):
#corrAI_SIDI <- cor.test(x=AI, y=SIDI, method = 'spearman')
#corrAI_SIDI

#Calculate the Spearman statistics (reported as mantel r) for all combinations of patch metrics
#mantel r = mantel coefficient
#pval1 = one-tailed p-value (null hypothesis: r <= 0).
#pval2 = one-tailed p-value (null hypothesis: r >= 0).
#pval3 = one-tailed p-value (null hypothesis: r = 0). #THIS IS THE ONE WE SHOULD CONSIDER
#llim = lower confidence limit
#ulim = upper confidence limit


corrAI_SIDI <- mantel(AI~SIDI, mrank=TRUE)
corrAI_SIEI <- mantel(AI~SIEI, mrank=TRUE)
corrAI_SHAPE_CV <- mantel(AI~SHAPE_CV, mrank=TRUE)
corrAI_FRAC_CV <- mantel(AI~FRAC_CV, mrank=TRUE)
corrAI_PARA_AM <- mantel(AI~PARA_AM, mrank=TRUE)
corrAI_CORE_CV <- mantel(AI~CORE_CV, mrank=TRUE)
corrAI_LPI <- mantel(AI~LPI, mrank=TRUE)
corrAI_DIVISION <- mantel(AI~DIVISION, mrank=TRUE)
corrAI_SPLIT <- mantel(AI~SPLIT, mrank=TRUE)

corrSIDI_SIEI <- mantel(SIDI~SIEI, mrank=TRUE)
corrSIDI_SHAPE_CV <- mantel(SIDI~SHAPE_CV, mrank=TRUE)
corrSIDI_FRAC_CV <- mantel(SIDI~FRAC_CV, mrank=TRUE)
corrSIDI_PARA_AM <- mantel(SIDI~PARA_AM, mrank=TRUE)
corrSIDI_CORE_CV <- mantel(SIDI~CORE_CV, mrank=TRUE)
corrSIDI_LPI <- mantel(SIDI~LPI, mrank=TRUE)
corrSIDI_DIVISION <- mantel(SIDI~DIVISION, mrank=TRUE)
corrSIDI_SPLIT <- mantel(SIDI~SPLIT, mrank=TRUE)

corrSIEI_SHAPE_CV <- mantel(SIEI~SHAPE_CV, mrank=TRUE)
corrSIEI_FRAC_CV <- mantel(SIEI~FRAC_CV, mrank=TRUE)
corrSIEI_PARA_AM <- mantel(SIEI~PARA_AM, mrank=TRUE)
corrSIEI_CORE_CV <- mantel(SIEI~CORE_CV, mrank=TRUE)
corrSIEI_LPI <- mantel(SIEI~LPI, mrank=TRUE)
corrSIEI_DIVISION <- mantel(SIEI~DIVISION, mrank=TRUE)
corrSIEI_SPLIT <- mantel(SIEI~SPLIT, mrank=TRUE)

corrSHAPE_CV_FRAC_CV <- mantel(SHAPE_CV~FRAC_CV, mrank=TRUE)
corrSHAPE_CV_PARA_AM <- mantel(SHAPE_CV~PARA_AM, mrank=TRUE)
corrSHAPE_CV_CORE_CV <- mantel(SHAPE_CV~CORE_CV, mrank=TRUE)
corrSHAPE_CV_LPI <- mantel(SHAPE_CV~LPI, mrank=TRUE)
corrSHAPE_CV_DIVISION <- mantel(SHAPE_CV~DIVISION, mrank=TRUE)
corrSHAPE_CV_SPLIT <- mantel(SHAPE_CV~SPLIT, mrank=TRUE)

corrFRAC_CV_PARA_AM <- mantel(FRAC_CV~PARA_AM, mrank=TRUE)
corrFRAC_CV_CORE_CV <- mantel(FRAC_CV~CORE_CV, mrank=TRUE)
corrFRAC_CV_LPI <- mantel(FRAC_CV~LPI, mrank=TRUE)
corrFRAC_CV_DIVISION <- mantel(FRAC_CV~DIVISION, mrank=TRUE)
corrFRAC_CV_SPLIT <- mantel(FRAC_CV~SPLIT, mrank=TRUE)

corrPARA_AM_CORE_CV <- mantel(PARA_AM~CORE_CV, mrank=TRUE)
corrPARA_AM_LPI <- mantel(PARA_AM~LPI, mrank=TRUE)
corrPARA_AM_DIVISION <- mantel(PARA_AM~DIVISION, mrank=TRUE)
corrPARA_AM_SPLIT <- mantel(PARA_AM~SPLIT, mrank=TRUE)

corrCORE_CV_LPI <- mantel(CORE_CV~LPI, mrank=TRUE)
corrCORE_CV_DIVISION <- mantel(CORE_CV~DIVISION, mrank=TRUE)
corrCORE_CV_SPLIT <- mantel(CORE_CV~SPLIT, mrank=TRUE)

corrLPI_DIVISION <- mantel(LPI~DIVISION, mrank=TRUE)
corrLPI_SPLIT <- mantel(LPI~SPLIT, mrank=TRUE)

corrDIVISION_SPLIT <- mantel(DIVISION~SPLIT, mrank=TRUE)


#print statistics

sink("SR_correlation_patchmetrix.txt")


cat("*****************************\n")
cat("SPEARMAN CORRELATION\n")
cat("*****************************\n")


cat("corrAI_SIDI\n")
cat("=============================\n")
corrAI_SIDI 
cat("=============================\n")
cat("\n")

cat("corrAI_SIEI\n")
cat("=============================\n")
corrAI_SIEI 
cat("=============================\n")
cat("\n")

cat("corrAI_SHAPE_CV\n")
cat("=============================\n")
corrAI_SHAPE_CV 
cat("=============================\n")
cat("\n")

cat("corrAI_FRAC_CV\n")
cat("=============================\n")
corrAI_FRAC_CV 
cat("=============================\n")
cat("\n")

cat("corrAI_PARA_AM\n")
cat("=============================\n")
corrAI_PARA_AM 
cat("=============================\n")
cat("\n")

cat("corrAI_CORE_CV\n")
cat("=============================\n")
corrAI_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrAI_LPI\n")
cat("=============================\n")
corrAI_LPI 
cat("=============================\n")
cat("\n")

cat("corrAI_DIVISION\n")
cat("=============================\n")
corrAI_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrAI_SPLIT\n")
cat("=============================\n")
corrAI_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SIEI\n")
cat("=============================\n")
corrSIDI_SIEI 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SHAPE_CV\n")
cat("=============================\n")
corrSIDI_SHAPE_CV 
cat("=============================\n")
cat("\n")

cat("corrSIDI_FRAC_CV\n")
cat("=============================\n")
corrSIDI_FRAC_CV 
cat("=============================\n")
cat("\n")

cat("corrSIDI_PARA_AM\n")
cat("=============================\n")
corrSIDI_PARA_AM 
cat("=============================\n")
cat("\n")

cat("corrSIDI_CORE_CV\n")
cat("=============================\n")
corrSIDI_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrSIDI_LPI\n")
cat("=============================\n")
corrSIDI_LPI 
cat("=============================\n")
cat("\n")

cat("corrSIDI_DIVISION\n")
cat("=============================\n")
corrSIDI_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SPLIT\n")
cat("=============================\n")
corrSIDI_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrSIEI_SHAPE_CV\n")
cat("=============================\n")
corrSIEI_SHAPE_CV 
cat("=============================\n")
cat("\n")

cat("corrSIEI_FRAC_CV\n")
cat("=============================\n")
corrSIEI_FRAC_CV 
cat("=============================\n")
cat("\n")

cat("corrSIEI_PARA_AM\n")
cat("=============================\n")
corrSIEI_PARA_AM 
cat("=============================\n")
cat("\n")

cat("corrSIEI_CORE_CV\n")
cat("=============================\n")
corrSIEI_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrSIEI_LPI\n")
cat("=============================\n")
corrSIEI_LPI 
cat("=============================\n")
cat("\n")

cat("corrSIEI_DIVISION\n")
cat("=============================\n")
corrSIEI_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrSIEI_SPLIT\n")
cat("=============================\n")
corrSIEI_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_FRAC_CV\n")
cat("=============================\n")
corrSHAPE_CV_FRAC_CV 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_PARA_AM\n")
cat("=============================\n")
corrSHAPE_CV_PARA_AM 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_CORE_CV\n")
cat("=============================\n")
corrSHAPE_CV_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_LPI\n")
cat("=============================\n")
corrSHAPE_CV_LPI 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_DIVISION\n")
cat("=============================\n")
corrSHAPE_CV_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_SPLIT\n")
cat("=============================\n")
corrSHAPE_CV_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_PARA_AM\n")
cat("=============================\n")
corrFRAC_CV_PARA_AM 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_CORE_CV\n")
cat("=============================\n")
corrFRAC_CV_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_LPI\n")
cat("=============================\n")
corrFRAC_CV_LPI 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_DIVISION\n")
cat("=============================\n")
corrFRAC_CV_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_SPLIT\n")
cat("=============================\n")
corrFRAC_CV_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_CORE_CV\n")
cat("=============================\n")
corrPARA_AM_CORE_CV 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_LPI\n")
cat("=============================\n")
corrPARA_AM_LPI 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_DIVISION\n")
cat("=============================\n")
corrPARA_AM_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_SPLIT\n")
cat("=============================\n")
corrPARA_AM_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_LPI\n")
cat("=============================\n")
corrCORE_CV_LPI 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_DIVISION\n")
cat("=============================\n")
corrCORE_CV_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_SPLIT\n")
cat("=============================\n")
corrCORE_CV_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrLPI_DIVISION\n")
cat("=============================\n")
corrLPI_DIVISION 
cat("=============================\n")
cat("\n")

cat("corrLPI_SPLIT\n")
cat("=============================\n")
corrLPI_SPLIT 
cat("=============================\n")
cat("\n")

cat("corrDIVISION_SPLIT\n")
cat("=============================\n")
corrDIVISION_SPLIT 
cat("=============================\n")
cat("\n")


#Calculate Pearson correlation between all patch metric combinations

corrAI_SIDI2 <- mantel(AI~SIDI, mrank=FALSE)
corrAI_SIEI2 <- mantel(AI~SIEI, mrank=FALSE)
corrAI_SHAPE_CV2 <- mantel(AI~SHAPE_CV, mrank=FALSE)
corrAI_FRAC_CV2 <- mantel(AI~FRAC_CV, mrank=FALSE)
corrAI_PARA_AM2 <- mantel(AI~PARA_AM, mrank=FALSE)
corrAI_CORE_CV2 <- mantel(AI~CORE_CV, mrank=FALSE)
corrAI_LPI2 <- mantel(AI~LPI, mrank=FALSE)
corrAI_DIVISION2 <- mantel(AI~DIVISION, mrank=FALSE)
corrAI_SPLIT2 <- mantel(AI~SPLIT, mrank=FALSE)

corrSIDI_SIEI2 <- mantel(SIDI~SIEI, mrank=FALSE)
corrSIDI_SHAPE_CV2 <- mantel(SIDI~SHAPE_CV, mrank=FALSE)
corrSIDI_FRAC_CV2 <- mantel(SIDI~FRAC_CV, mrank=FALSE)
corrSIDI_PARA_AM2 <- mantel(SIDI~PARA_AM, mrank=FALSE)
corrSIDI_CORE_CV2 <- mantel(SIDI~CORE_CV, mrank=FALSE)
corrSIDI_LPI2 <- mantel(SIDI~LPI, mrank=FALSE)
corrSIDI_DIVISION2 <- mantel(SIDI~DIVISION, mrank=FALSE)
corrSIDI_SPLIT2 <- mantel(SIDI~SPLIT, mrank=FALSE)

corrSIEI_SHAPE_CV2 <- mantel(SIEI~SHAPE_CV, mrank=FALSE)
corrSIEI_FRAC_CV2 <- mantel(SIEI~FRAC_CV, mrank=FALSE)
corrSIEI_PARA_AM2 <- mantel(SIEI~PARA_AM, mrank=FALSE)
corrSIEI_CORE_CV2 <- mantel(SIEI~CORE_CV, mrank=FALSE)
corrSIEI_LPI2 <- mantel(SIEI~LPI, mrank=FALSE)
corrSIEI_DIVISION2 <- mantel(SIEI~DIVISION, mrank=FALSE)
corrSIEI_SPLIT2 <- mantel(SIEI~SPLIT, mrank=FALSE)

corrSHAPE_CV_FRAC_CV2 <- mantel(SHAPE_CV~FRAC_CV, mrank=FALSE)
corrSHAPE_CV_PARA_AM2 <- mantel(SHAPE_CV~PARA_AM, mrank=FALSE)
corrSHAPE_CV_CORE_CV2 <- mantel(SHAPE_CV~CORE_CV, mrank=FALSE)
corrSHAPE_CV_LPI2 <- mantel(SHAPE_CV~LPI, mrank=FALSE)
corrSHAPE_CV_DIVISION2 <- mantel(SHAPE_CV~DIVISION, mrank=FALSE)
corrSHAPE_CV_SPLIT2 <- mantel(SHAPE_CV~SPLIT, mrank=FALSE)

corrFRAC_CV_PARA_AM2 <- mantel(FRAC_CV~PARA_AM, mrank=FALSE)
corrFRAC_CV_CORE_CV2 <- mantel(FRAC_CV~CORE_CV, mrank=FALSE)
corrFRAC_CV_LPI2 <- mantel(FRAC_CV~LPI, mrank=FALSE)
corrFRAC_CV_DIVISION2 <- mantel(FRAC_CV~DIVISION, mrank=FALSE)
corrFRAC_CV_SPLIT2 <- mantel(FRAC_CV~SPLIT, mrank=FALSE)

corrPARA_AM_CORE_CV2 <- mantel(PARA_AM~CORE_CV, mrank=FALSE)
corrPARA_AM_LPI2 <- mantel(PARA_AM~LPI, mrank=FALSE)
corrPARA_AM_DIVISION2 <- mantel(PARA_AM~DIVISION, mrank=FALSE)
corrPARA_AM_SPLIT2 <- mantel(PARA_AM~SPLIT, mrank=FALSE)

corrCORE_CV_LPI2 <- mantel(CORE_CV~LPI, mrank=FALSE)
corrCORE_CV_DIVISION2 <- mantel(CORE_CV~DIVISION, mrank=FALSE)
corrCORE_CV_SPLIT2 <- mantel(CORE_CV~SPLIT, mrank=FALSE)

corrLPI_DIVISION2 <- mantel(LPI~DIVISION, mrank=FALSE)
corrLPI_SPLIT2 <- mantel(LPI~SPLIT, mrank=FALSE)

corrDIVISION_SPLIT2 <- mantel(DIVISION~SPLIT, mrank=FALSE)


#print statistics


cat("*****************************\n")
cat("PEARSONS CORRELATION\n")
cat("*****************************\n")

cat("corrAI_SIDI\n")
cat("=============================\n")
corrAI_SIDI2 
cat("=============================\n")
cat("\n")

cat("corrAI_SIEI\n")
cat("=============================\n")
corrAI_SIEI2 
cat("=============================\n")
cat("\n")

cat("corrAI_SHAPE_CV\n")
cat("=============================\n")
corrAI_SHAPE_CV2 
cat("=============================\n")
cat("\n")

cat("corrAI_FRAC_CV\n")
cat("=============================\n")
corrAI_FRAC_CV2 
cat("=============================\n")
cat("\n")

cat("corrAI_PARA_AM\n")
cat("=============================\n")
corrAI_PARA_AM2 
cat("=============================\n")
cat("\n")

cat("corrAI_CORE_CV\n")
cat("=============================\n")
corrAI_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrAI_LPI\n")
cat("=============================\n")
corrAI_LPI2 
cat("=============================\n")
cat("\n")

cat("corrAI_DIVISION\n")
cat("=============================\n")
corrAI_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrAI_SPLIT\n")
cat("=============================\n")
corrAI_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SIEI\n")
cat("=============================\n")
corrSIDI_SIEI2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SHAPE_CV\n")
cat("=============================\n")
corrSIDI_SHAPE_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_FRAC_CV\n")
cat("=============================\n")
corrSIDI_FRAC_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_PARA_AM\n")
cat("=============================\n")
corrSIDI_PARA_AM2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_CORE_CV\n")
cat("=============================\n")
corrSIDI_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_LPI\n")
cat("=============================\n")
corrSIDI_LPI2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_DIVISION\n")
cat("=============================\n")
corrSIDI_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrSIDI_SPLIT\n")
cat("=============================\n")
corrSIDI_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_SHAPE_CV\n")
cat("=============================\n")
corrSIEI_SHAPE_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_FRAC_CV\n")
cat("=============================\n")
corrSIEI_FRAC_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_PARA_AM\n")
cat("=============================\n")
corrSIEI_PARA_AM2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_CORE_CV\n")
cat("=============================\n")
corrSIEI_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_LPI\n")
cat("=============================\n")
corrSIEI_LPI2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_DIVISION\n")
cat("=============================\n")
corrSIEI_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrSIEI_SPLIT\n")
cat("=============================\n")
corrSIEI_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_FRAC_CV\n")
cat("=============================\n")
corrSHAPE_CV_FRAC_CV2 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_PARA_AM\n")
cat("=============================\n")
corrSHAPE_CV_PARA_AM2
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_CORE_CV\n")
cat("=============================\n")
corrSHAPE_CV_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_LPI\n")
cat("=============================\n")
corrSHAPE_CV_LPI2
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_DIVISION\n")
cat("=============================\n")
corrSHAPE_CV_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrSHAPE_CV_SPLIT\n")
cat("=============================\n")
corrSHAPE_CV_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_PARA_AM\n")
cat("=============================\n")
corrFRAC_CV_PARA_AM2 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_CORE_CV\n")
cat("=============================\n")
corrFRAC_CV_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_LPI\n")
cat("=============================\n")
corrFRAC_CV_LPI2 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_DIVISION\n")
cat("=============================\n")
corrFRAC_CV_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrFRAC_CV_SPLIT\n")
cat("=============================\n")
corrFRAC_CV_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_CORE_CV\n")
cat("=============================\n")
corrPARA_AM_CORE_CV2 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_LPI\n")
cat("=============================\n")
corrPARA_AM_LPI2 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_DIVISION\n")
cat("=============================\n")
corrPARA_AM_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrPARA_AM_SPLIT\n")
cat("=============================\n")
corrPARA_AM_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_LPI\n")
cat("=============================\n")
corrCORE_CV_LPI2 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_DIVISION\n")
cat("=============================\n")
corrCORE_CV_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrCORE_CV_SPLIT\n")
cat("=============================\n")
corrCORE_CV_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrLPI_DIVISION\n")
cat("=============================\n")
corrLPI_DIVISION2 
cat("=============================\n")
cat("\n")

cat("corrLPI_SPLIT\n")
cat("=============================\n")
corrLPI_SPLIT2 
cat("=============================\n")
cat("\n")

cat("corrDIVISION_SPLIT\n")
cat("=============================\n")
corrDIVISION_SPLIT2 
cat("=============================\n")
cat("\n")

sink()

############################################################################################################################
## Richard Working code                                            #########################################################
############################################################################################################################
df1<-combined_dataset#I'm lazy...
View(df1)
##For my poor little brain: 
#lab 12 has ID which corresponds to landscape, Name:NA?,LogDc.km:Rst,LENGTH:euclid_dist,next several are the variables
#Site1=[pop1 and site 2=pop2
##Our Rst values are non-normal...does that matter? It might. I did it for fun:
df1$trans.RST<-log(df1$RST)#introduces NaNs, not using it.
Zl <- lapply(c("site_1","site_2"), function(nm) 
  Matrix:::fac2sparse(df1[[nm]], "d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])

##############################################################
#using continuous metrics first
############################################################

mod1 <- lme4::lFormula(RST ~ Sa+S10z+Ssk+Sku+Sdr+Sbi+Std+Stdi+Sfd+Srwi+(1|site_1), data = df1, REML = TRUE)
#Like Lab 12 got some scale issues, check where:
head(df1[,3:12])#lookslike s10z std and sdr
S10zz<-scale(df1$S10z,center=T,scale=T)
Sdrz<-scale(df1$Sdr,center=T,scale=T)
Stdz<-scale(df1$Std,center=T,scale=T)
df1<-cbind(df1,S10zz)
df1<-cbind(df1,Sdrz)
df1<-cbind(df1,Stdz)
head(df1[,c(3:12,36,37)])
##Cool...retry model (switching z scaled variables in) before I mess more.
mod1 <- lme4::lFormula(RST ~ Sa+S10zz+Ssk+Sku+Sdrz+Sbi+Stdz+Stdi+Sfd+Srwi+(1|site_1), data = df1, REML = TRUE)
##still funky, checking:
head(df1[,c(3,5:6,8,10:12,36,37,38,39)])
##scaling Sa
df1$Saz<-scale(df1$Sa,center=T,scale=T)
#Now we're using these:
head(df1[,c(5:6,8,10:12,36,37,38,39)])#and if Sfd causes problems, so help me God, I'll delete it...
mod1 <- lme4::lFormula(RST ~ Saz+S10zz+Ssk+Sku+Sdrz+Sbi+Stdz+Stdi+Sfd+Srwi+(1|site_1), data = df1, REML = TRUE)
##Still having some issues, but going to tolerate them for now.Testing multicollinearity
df1.prevars <- with(df1, data.frame(Saz,S10zz,Ssk,Sku,Sdrz,Sbi,Stdz,Stdi,Sfd,Srwi))
usdm::vif(df1.prevars)
#Great Scot: Skuis pulled first
df1.vars <- with(df1, data.frame(Saz,S10zz,Ssk,Sdrz,Sbi,Stdz,Stdi,Sfd,Srwi))
usdm::vif(df1.vars)
#Much happier, all less than 10. I'm considering that acceptable.
#let's update the model
Zl <- lapply(c("site_1","site_2"), function(nm) 
  Matrix:::fac2sparse(CSFdata[[nm]], "d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])
mod1 <- lme4::lFormula(RST ~ Saz+S10zz+Ssk+Sdrz+Sbi+Stdz+Stdi+Sfd+Srwi+(1|site_1), data = df1, REML = TRUE)
##woot woot! no scale issues

#################################################PROBLEM WITH THIS LINE################################
mod1$reTrms$Zt <- ZZ
#############################################################################
names(df1)
##refit model
dfun <- do.call(lme4::mkLmerDevfun, mod1)
opt <- lme4::optimizeLmer(dfun)
mod1 <- lme4::mkMerMod(environment(dfun), opt, mod1$reTrms, fr = mod1$fr)
summary(mod1)
#Plot it
plot(mod1)
#We have a trumpet shape...large variance at large numbers.Multiplicative errors with increasing dist?
hist(residuals(mod1))


##MLPE FUNC
MLPE <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = TRUE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

##Defining variables
mod2 <- MLPE(RST ~ Ssk + Sbi + Stdi + Sfd+Srwi+S10zz+Sdrz+Stdz+Saz + (1|site_1), df1)
mod3 <- MLPE(RST ~ Ssk + Sbi + Sfd + S10zz +  Sdrz + Stdz + Saz + (1|site_1), df1)
mod4 <- MLPE(RST ~ Ssk + Sfd + S10zz + Sdrz + Stdz + (1|site_1), df1)

MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

mod1noREML <- MLPEnoREML (RST ~ Saz+S10zz+Ssk+Sdrz+Sbi+Stdz+Stdi+Sfd+Srwi+(1|site_1), df1)
mod2noREML <- MLPEnoREML(RST ~ Ssk + Sbi + Stdi + Sfd+Srwi+S10zz+Sdrz+Stdz+Saz + (1|site_1), df1)
mod3noREML <- MLPEnoREML(RST ~ Ssk + Sbi + Sfd + S10zz +  Sdrz + Stdz + Saz + (1|site_1), df1)
mod4noREML <- MLPEnoREML(RST ~ Ssk + Sfd + S10zz + Sdrz + Stdz + (1|site_1), df1)

Models <- list(Full=mod1noREML, fullagain=mod2noREML, thinned=mod3noREML, 
               reduced=mod4noREML)
CSF.IC <- data.frame(AIC = sapply(Models, AIC),
                     BIC = sapply(Models, BIC)) 
CSF.IC <- data.frame(CSF.IC, k = sapply(Models, function(ls) attr(logLik(ls), "df")))
CSF.IC
CSF.IC$AICc <- CSF.IC$AIC + 2*CSF.IC$k*(CSF.IC$k+1)/(48-CSF.IC$k-1)
CSF.IC##has corrected
##Model weights
AICcmin <- min(CSF.IC$AICc)
RL <- exp(-0.5*(CSF.IC$AICc - AICcmin))
sumRL <- sum(RL)
CSF.IC$AICcmin <- RL/sumRL
BICmin <- min(CSF.IC$BIC)
RL.B <- exp(-0.5*(CSF.IC$BIC - BICmin))
sumRL.B <- sum(RL.B)
CSF.IC$BICew <- RL.B/sumRL.B
round(CSF.IC,3)

##CIs eventually:
ModelsREML <- list(Full=mod1, full2=mod2, thinned=mod3,
                   reduced=mod4)
confint(ModelsREML$thinned, level = 0.95, method = "Wald")
confint(ModelsREML$Full, level = 0.95, method = "Wald")
confint(ModelsREML$reduced, level = 0.95, method = "Wald")


################################################################################
#Patch metrics
##############################################################################
modA <- lme4::lFormula(RST ~ AI+SIDI+SIEI+AREA_CV+SHAPE_CV+FRAC_CV+PARA_AM+CORE_CV+LPI+DIVISION+SPLIT+ED+LSI+(1|site_1), data = df1, REML = TRUE)
#Like Lab 12 got some scale issues, check where:
head(df1[,19:31])
df1$AIz<-scale(df1$AI,center=T,scale=T)
df1$AREA_CVz<-scale(df1$AREA_CV,center=T,scale=T)
df1$SHAPEz<-scale(df1$SHAPE_CV,center=T,scale=T)
df1$CORE_CVz<-scale(df1$CORE_CV,center=T,scale=T)
df1$LPIz<-scale(df1$LPI,center=T,scale=T)
df1$LSIz<-scale(df1$LSI,center=T,scale=T)
head(df1[,c(19:31,40:45)])
##Cool...retry model (switching z scaled variables in) before I mess more.
modA <- lme4::lFormula(RST ~ AIz+SIDI+SIEI+AREA_CVz+SHAPEz+FRAC_CV+PARA_AM+CORE_CVz+LPIz+DIVISION+SPLIT+ED+LSIz+(1|site_1), data = df1, REML = TRUE)

##Testing multicollinearity

df1.patchvars <- with(df1, data.frame(AIz,SIDI,SIEI,AREA_CVz,SHAPEz,FRAC_CV,PARA_AM,CORE_CVz,LPIz,DIVISION,SPLIT,ED,LSIz))
usdm::vif(df1.patchvars)

#Great Scot: They're astronomical.

#I reduced them one variable at a time until I got below ten.  These are the lucky winners.
df1.patchvars <- with(df1, data.frame(FRAC_CV,SPLIT,ED,LSIz))
usdm::vif(df1.patchvars)

Zl <- lapply(c("site_1","site_2"), function(nm) 
  Matrix:::fac2sparse(CSFdata[[nm]], "d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])
modA <- lme4::lFormula(RST ~ FRAC_CV+SPLIT+ED+LSIz+(1|site_1), data = df1, REML = TRUE)


#################################################PROBLEM WITH THIS LINE################################
mod1$reTrms$Zt <- ZZ
#############################################################################

##refit model
dfun <- do.call(lme4::mkLmerDevfun, modA)
opt <- lme4::optimizeLmer(dfun)
modA <- lme4::mkMerMod(environment(dfun), opt, modA$reTrms, fr = modA$fr)
summary(modA)
#Plot it
plot(modA)
#We have a trumpet shape...large variance at large numbers.Multiplicative errors with increasing dist?
hist(residuals(modA))


##Defining variables
modB <- MLPE(RST ~ SPLIT+ED+LSIz + (1|site_1), df1)
modC <- MLPE(RST ~ SPLIT+ED + (1|site_1), df1)

modAnoREML <- MLPEnoREML (RST ~ FRAC_CV+SPLIT+ED+LSIz+(1|site_1), df1)
modBnoREML <- MLPEnoREML(RST ~ SPLIT+ED+LSIz + (1|site_1), df1)
modCnoREML <- MLPEnoREML(RST ~ SPLIT+ED + (1|site_1), df1)

Models <- list(Full=modAnoREML, stepwise1=modBnoREML, stepwise3=modCnoREML)
CSF.IC <- data.frame(AIC = sapply(Models, AIC),
                     BIC = sapply(Models, BIC)) 
CSF.IC <- data.frame(CSF.IC, k = sapply(Models, function(ls) attr(logLik(ls), "df")))
CSF.IC
CSF.IC$AICc <- CSF.IC$AIC + 2*CSF.IC$k*(CSF.IC$k+1)/(48-CSF.IC$k-1)
CSF.IC##has corrected
##Model weights
AICcmin <- min(CSF.IC$AICc)
RL <- exp(-0.5*(CSF.IC$AICc - AICcmin))
sumRL <- sum(RL)
CSF.IC$AICcmin <- RL/sumRL
BICmin <- min(CSF.IC$BIC)
RL.B <- exp(-0.5*(CSF.IC$BIC - BICmin))
sumRL.B <- sum(RL.B)
CSF.IC$BICew <- RL.B/sumRL.B
round(CSF.IC,3)

##CIs eventually:
ModelsREML <- list(Full=modA, SW1=modB, SW2=modC)
confint(ModelsREML$Full, level = 0.95, method = "Wald")
confint(ModelsREML$SW1, level = 0.95, method = "Wald")
confint(ModelsREML$SW2, level = 0.95, method = "Wald")