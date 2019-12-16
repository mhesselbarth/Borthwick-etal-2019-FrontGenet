# Testing best models against Rst values to calculate kappa statistics #

# before this script runs, you need to run 05_run_models.R and 06_LCP_models.R to 
# have all required data in workspace

# load libraries
library(irr)
library(rhoR)
library(psych)
library(rminer)

# key out models starting with lsm:
p <- predict(model_fit_lsm, df_lsm)
a <- df_lsm$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
lsm_rho <- rho(abs(k[[2]]),OcSBaserate = b[[2]],
               testSetLength = nrow(d), ScSKappaThreshold = 0,ScSKappaMin = 0)

#then for surface:
p <- predict(model_fit_surface, df_surface)
a <- df_surface$DPS
d <- cbind(a,p)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
surf_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]],
                testSetLength = nrow(d), ScSKappaThreshold = 0, ScSKappaMin = 0)

#lsm with NA:
p <- predict(model_fit_lsm_NA, df_lsm_NA)
a <- df_lsm_NA$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
NA_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d),
              ScSKappaThreshold = 0, ScSKappaMin = 0)

#LCP
p <- predict(model_fit_lcp, df_lcp)
a <- df_lcp$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
lcp_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d), 
               ScSKappaThreshold = 0, ScSKappaMin = 0)

#RES
p <- predict(model_fit_res, df_res)
a <- df_res$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
res.rho <- rho(abs(k[[2]]),OcSBaserate = b[[2]], testSetLength = nrow(d), 
               ScSKappaThreshold = 0,ScSKappaMin = 0)


#IBD
p <- predict(model_fit_ibd, df_ibd)
a <- df_ibd$DPS
d <- cbind(a,p)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
ibd_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d), 
               ScSKappaThreshold = 0,ScSKappaMin = 0)
