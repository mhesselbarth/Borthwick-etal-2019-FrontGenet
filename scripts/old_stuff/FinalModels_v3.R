############################################################################################################################
## Richard Working code                                            #########################################################
############################################################################################################################
library(ecodist)
require(ecodist)
library(readr)
require(readr)
library(LandGenCourse)
library(ggplot2)
library(MuMIn)
combined_dataset <- read_csv("combined_datasetV3.csv") 

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

mod1 <- lme4::lFormula(RST ~ Sa+S10z+Ssk+Sku+Sdr+Sbi+Std+Stdi+Sfd+Srwi+Geo_dist+(1|site_1), data = df1, REML = TRUE)
#Like Lab 12 got some scale issues, check where:
head(df1[,3:12])#lookslike s10z std and sdr
##Standardize EVERYTHING...
df1z<-df1
head(df1z)
df1z[,c(3:14,17:31)]<-scale(df1z[,c(3:14,17:31)],center = T,scale=T)
##Cool...retry model (switching z scaled variables in) before I mess more.
mod1 <- lme4::lFormula(RST ~ Sa+S10z+Ssk+Sku+Sdr+Sbi+Std+Stdi+Sfd+Srwi+Geo_dist+(1|site_1), data = df1z, REML = TRUE)

##Testing multicollinearity
df1.prevars <- with(df1z, data.frame(Sa,S10z,Ssk,Sku,Sdr,Sbi,Std,Stdi,Sfd,Srwi))
usdm::vif(df1.prevars)
#Great Scot: Sku is pulled first
df1.vars <- with(df1z, data.frame(Sa,S10z,Ssk,Sdr,Sbi,Std,Stdi,Sfd,Srwi))
usdm::vif(df1.vars)
#Much happier, all less than 10. I'm considering that acceptable.
#let's update the model
Zl <- lapply(c("site_1","site_2"), function(nm) 
  Matrix:::fac2sparse(CSFdata[[nm]], "d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])
mod1 <- lme4::lFormula(RST ~ Sa+S10z+Ssk+Sdr+Sbi+Std+Stdi+Sfd+Srwi+(1|site_1), data = df1z, REML = TRUE)


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
mod2 <- MLPE(RST ~ Ssk + Sbi + Stdi + Sfd+Srwi+S10z+Sdr+Std+Sa + (1|site_1), df1z)
mod3 <- MLPE(RST ~ Ssk + Sbi + Sfd + S10z +  Sdr + Std + Sa + (1|site_1), df1z)
mod4 <- MLPE(RST ~ Ssk + Sfd + S10z + Sdr + Std + (1|site_1), df1z)

MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  ##mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

mod1noREML <- MLPEnoREML (RST ~ Sa+S10z+Ssk+Sdr+Sbi+Std+Stdi+Sfd+Srwi+(1|site_1), df1z)
mod2noREML <- MLPEnoREML(RST ~ Ssk + Sbi + Stdi + Sfd+Srwi+S10z+Sdr+Std+Sa + (1|site_1), df1z)
mod3noREML <- MLPEnoREML(RST ~ Ssk + Sbi + Sfd + S10z +  Sdr + Stdz + Sa + (1|site_1), df1z)
mod4noREML <- MLPEnoREML(RST ~ Ssk + Sfd + S10z + Sdr + Std + (1|site_1), df1z)

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
modA <- lme4::lFormula(RST ~ AI+SIDI+SIEI+AREA_CV+SHAPE_CV+FRAC_CV+PARA_AM+CORE_CV+LPI+DIVISION+SPLIT+ED+LSI+(1|site_1), data = df1z, REML = TRUE)
#Like Lab 12 got some scale issues, check where:
modA <- lme4::lFormula(RST ~ AI+SIDI+SIEI+AREA_CV+SHAPE_CV+FRAC_CV+PARA_AM+CORE_CV+LPI+DIVISION+SPLIT+ED+LSI+(1|site_1), data = df1z, REML = TRUE)

##Testing multicollinearity

df1.patchvars <- with(df1z, data.frame(AI,SIDI,SIEI,AREA_CV,SHAPE_CV,FRAC_CV,PARA_AM,CORE_CV,LPI,DIVISION,SPLIT,ED,LSI))
usdm::vif(df1.patchvars)

#Great Scot: They're astronomical.

#I reduced them one variable at a time until I got below ten.  These are the lucky winners.
df1.patchvars <- with(df1z, data.frame(FRAC_CV,SPLIT,ED,LSI))
usdm::vif(df1.patchvars)

Zl <- lapply(c("site_1","site_2"), function(nm) 
  Matrix:::fac2sparse(CSFdata[[nm]], "d", drop=FALSE))
ZZ <- Reduce("+", Zl[-1], Zl[[1]])
modA <- lme4::lFormula(RST ~ FRAC_CV+SPLIT+ED+LSI+(1|site_1), data = df1, REML = TRUE)


#################################################PROBLEM WITH THIS LINE################################
##mod1$reTrms$Zt <- ZZ
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
modB <- MLPE(RST ~ SPLIT+ED+LSI + (1|site_1), df1z)
modC <- MLPE(RST ~ LSI + (1|site_1), df1z)

modAnoREML <- MLPEnoREML (RST ~ FRAC_CV+SPLIT+ED+LSI+(1|site_1), df1z)
modBnoREML <- MLPEnoREML(RST ~ SPLIT+ED+LSI + (1|site_1), df1z)
modCnoREML <- MLPEnoREML(RST ~ LSI + (1|site_1), df1z)

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

####INSTEAD OF THE ABOVE USE DREDGE#####
# options('na.actions'="na.omit")
pmm.model <- lme4::lmer(RST ~ FRAC_CV+SPLIT+ED+LSI+Geo_dist+(1|site_1), 
                  data = df1z, REML=F, 
                  na.action = 'na.fail')
output <- dredge(pmm.model)
output ####Fixed term is "(Intercept)"###

# options('na.actions'="na.omit")
sm.model <- lme4::lmer(RST ~ Sa+S10z+Ssk+Sku+Sdr+Sbi+Std+Stdi+Sfd+Srwi+Geo_dist+(1|site_1), 
                 data = df1z, REML = F, 
                 na.action = 'na.fail')
output2 <- dredge(sm.model)
output2 ####Fixed term is "(Intercept)"###
