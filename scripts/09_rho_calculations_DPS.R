# Testing best models against Rst values to calculate kappa statistics #

# before this script runs, you need to run 05_run_models.R and 06_LCP_models.R to 
# have all required data in workspace

# load libraries
library(irr)
library(rhoR)
library(psych)
library(rminer)

# key out models starting with lsm:
model_fit_lsm <- modular_function(response = "DPS", 
                                  explanatory = best_model_lsm[[1]],
                                  random =  "(1|site_a)", 
                                  data = df_lsm, 
                                  REML = TRUE, 
                                  ZZ = ZZ_landscape,
                                  verbose = FALSE)

p <- predict(model_fit_lsm, df_lsm)
a <- df_lsm$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
lsm_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d), 
               ScSKappaThreshold = k$confid[2,3], ScSKappaMin = k$confid[2,1])

#then for surface:
model_fit_surface <- modular_function(response = "DPS", 
                                      explanatory = best_model_surface[[1]],
                                      random =  "(1|site_1)", 
                                      data = df_surface, 
                                      REML = TRUE, 
                                      ZZ = ZZ_landscape,
                                      verbose = FALSE)

p <- predict(model_fit_surface, df_surface)
a <- df_surface$DPS
d <- cbind(a,p)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
surf_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d),
              ScSKappaThreshold = k$confid[2,3], ScSKappaMin = k$confid[2,1])

#lsm with NA:
model_fit_lsm_NA <- modular_function(response = "DPS", 
                                     explanatory = best_model_lsm_NA[[1]],
                                     random =  "(1|site_a)", 
                                     data = df_lsm_NA, 
                                     REML = TRUE, 
                                     ZZ = ZZ_landscape,
                                     verbose = FALSE)

p <- predict(model_fit_lsm_NA, df_lsm_NA)
a <- df_lsm_NA$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
NA_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d),
              ScSKappaThreshold = k$confid[2,3], ScSKappaMin = k$confid[2,1])

#LCP
model_fit_lcp <- modular_function(response = "DPS", 
                                  explanatory = "least_cost_scaled",
                                  random =  "(1|site_1)", 
                                  data = df_lcp, 
                                  REML = TRUE, 
                                  ZZ = ZZ_landscape,
                                  verbose = FALSE)

p <- predict(model_fit_lcp, df_lcp)
a <- df_lcp$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
lcp_rho <- rho(abs(k[[2]]),OcSBaserate = b[[2]], testSetLength = nrow(d),
               ScSKappaThreshold = k$confid[2,3],ScSKappaMin = k$confid[2,1])

#RES
model_fit_res <- modular_function(response = "DPS", 
                                  explanatory = "resistance_scaled",
                                  random =  "(1|site_1)", 
                                  data = df_res, 
                                  REML = TRUE, 
                                  ZZ = ZZ_landscape,
                                  verbose = FALSE)

p <- predict(model_fit_res, df_res)
a <- df_res$DPS
d <- cbind(p,a)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
res.rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d),
               ScSKappaThreshold = k$confid[2,3], ScSKappaMin = k$confid[2,1])


#IBD
model_fit_ibd <- modular_function(response = "DPS", 
                                  explanatory = "dist_scaled",
                                  random =  "(1|site_1)", 
                                  data = df_ibd, 
                                  REML = TRUE, 
                                  ZZ = ZZ_landscape,
                                  verbose = FALSE)

p <- predict(model_fit_ibd, df_ibd)
a <- df_ibd$DPS
d <- cbind(a,p)
k <- cohen.kappa(d)
b <- baserate(d = contingencyTable)
ibd_rho <- rho(abs(k[[2]]), OcSBaserate = b[[2]], testSetLength = nrow(d),
               ScSKappaThreshold = k$confid[2,3], ScSKappaMin = k$confid[2,1])

# results
tibble::tibble(
  surface = surf_rho,
  lsm = lsm_rho,
  lsm_na = NA_rho,
  resistance = res.rho,
  ibd = ibd_rho)
