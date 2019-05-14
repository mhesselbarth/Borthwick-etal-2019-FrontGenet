# load libraries
library(tidyverse)
library(lme4)
library(Matrix)
library(MuMIn)
library(usdn)

# model optimization function
modular_function <- function(variables, data, REML = TRUE, ZZ) {

  # parse the data and formula
  model <- lme4::lFormula(variables, data = data, REML = REML)
  
  # replace ZZ matrix
  model$reTrms$Zt <- ZZ

  # create the deviance function to be optimized
  deviance_function <- do.call(lme4::mkLmerDevfun, model)

  # optimize the deviance function:
  optimization <- lme4::optimizeLmer(deviance_function)

  # package up the results
  model <- lme4::mkMerMod(rho = environment(deviance_function),
                          opt = optimization,
                          reTrms = model$reTrms,
                          fr = model$fr)
  
  return(model)
}

#### Surface metrics ####

# import data surface metrics
surface_metrics <- readr::read_rds("data/Output/surface_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add RST to data sets
surface_metrics <- dplyr::left_join(surface_metrics, rst, by = c("site_1", "site_2"))

# Create Zl and ZZ matrix
Zl_surface <- lapply(c("site_1","site_2"), function(x) {Matrix::fac2sparse(surface_metrics[[x]], "d", drop = FALSE)})

ZZ_surface <- Reduce("+", Zl_surface[-1], Zl_surface[[1]])

# model surface metrics
surface_metrics_model <- lme4::lFormula(RST ~ Sa + S10z + Ssk + Sku + Sdr + Sbi + 
                                          Std + Stdi + Sfd + Srwi + (1|site_1), 
                                        data = surface_metrics, REML = TRUE)

# Warning message: Some predictor variables are on 
# very different scales: consider rescaling 

# rescale some metrics
surface_metrics <- dplyr::mutate(surface_metrics, 
                                 S10z_scaled = as.numeric(scale(S10z, center = TRUE, scale = TRUE)),
                                 Sdr_scaled = as.numeric(scale(Sdr, center = TRUE, scale = TRUE)), 
                                 Std_scaled = as.numeric(scale(Std, center = TRUE, scale = TRUE)),
                                 Sa_scaled = as.numeric(scale(Sa, center = TRUE, scale = TRUE)))

# rerun model
surface_metrics_model <- lme4::lFormula(RST ~ Sa_scaled + S10z_scaled + Ssk + Sku + Sdr_scaled + Sbi + 
                                          Std_scaled + Stdi + Sfd + Srwi + (1|site_1), 
                                        data = surface_metrics, REML = TRUE)

# Checking for variance inflation multicollinearity
dplyr::select(surface_metrics, 
              Sa_scaled, S10z_scaled, Ssk, Sku, Sdr_scaled, Sbi, Std_scaled, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vif()

# removing Sku
dplyr::select(surface_metrics, 
              Sa_scaled, S10z_scaled, Ssk, Sdr_scaled, Sbi, Std_scaled, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vif()

# update model
surface_metrics_model_1 <- modular_function(RST ~ Sa_scaled + S10z_scaled + Ssk + 
                                              Sdr_scaled + Sbi + Std_scaled + 
                                              Stdi + Sfd + Srwi + (1|site_1), 
                                            data = surface_metrics, 
                                            ZZ = ZZ_surface)

# look at model summary
summary(surface_metrics_model_1)

# plot model results
ggplot() + 
  geom_point(aes(x = fitted(surface_metrics_model_1), 
                 y = residuals(surface_metrics_model_1)), 
             pch = 1, size = 2.5) + 
  geom_hline(yintercept = 0) + 
  labs(x = "fitted values", y = "residuals") +
  theme_bw()

# Defining variables for different models
surface_metrics_model_2 <- modular_function(RST ~ Ssk + Sbi + Stdi + Sfd + Srwi + 
                                              S10z_scaled + Sdr_scaled + Std_scaled + 
                                              Sa_scaled + (1|site_1), 
                                            data = surface_metrics,
                                            ZZ = ZZ_surface)

surface_metrics_model_3 <- modular_function(RST ~ Ssk + Sbi + Sfd + S10z_scaled + 
                                              Sdr_scaled + Std_scaled + 
                                              Sa_scaled + (1|site_1), 
                                            data = surface_metrics,
                                            ZZ = ZZ_surface)

surface_metrics_model_4 <- modular_function(RST ~ Ssk + Sfd + S10z_scaled + 
                                              Sdr_scaled + Std_scaled + (1|site_1),
                                            data = surface_metrics,
                                            ZZ = ZZ_surface)

# setting REML to FALSE
surface_metrics_model_1_nr <- modular_function(RST ~ Sa_scaled + S10z_scaled + Ssk + 
                                                 Sdr_scaled + Sbi + Std_scaled + 
                                                 Stdi + Sfd + Srwi + (1|site_1), 
                                               data = surface_metrics, REML = FALSE,
                                               ZZ = ZZ_surface)

surface_metrics_model_2_nr <- modular_function(RST ~ Ssk + Sbi + Stdi + Sfd + Srwi + 
                                                 S10z_scaled + Sdr_scaled + Std_scaled + 
                                                 Sa_scaled + (1|site_1), 
                                               data = surface_metrics, REML = FALSE,
                                               ZZ = ZZ_surface)

surface_metrics_model_3_nr <- modular_function(RST ~ Ssk + Sbi + Sfd + S10z_scaled + 
                                                 Sdr_scaled + Std_scaled + 
                                                 Sa_scaled + (1|site_1), 
                                               data = surface_metrics, REML = FALSE,
                                               ZZ = ZZ_surface)

surface_metrics_model_4_nr <- modular_function(RST ~ Ssk + Sfd + S10z_scaled + 
                                                 Sdr_scaled + Std_scaled + (1|site_1),
                                               data = surface_metrics, REML = FALSE,
                                               ZZ = ZZ_surface)

# put all models in list
models_list_surface_nr <- list(model_1 = surface_metrics_model_1_nr, 
                               model_2 = surface_metrics_model_2_nr, 
                               model_3 = surface_metrics_model_3_nr, 
                               model_4 = surface_metrics_model_4_nr)

# calculating AIC and BIC
information_criterion_surface_nr <- purrr::map_dfr(models_list_surface_nr, function(x) {
  data.frame(AIC = AIC(x), 
             BIC = BIC(x),
             log_like = attr(logLik(x), "df"))
  }, .id = "model")

# calculating corrected AICc and BICc
information_criterion_surface_nr <- dplyr::mutate(information_criterion_surface_nr, 
                                                  AICc = AIC + 2 * log_like *
                                                    (log_like + 1) / (48 - log_like - 1),
                                                  AICc_min = exp(-0.5 * (AICc - min(AICc))) / sum(exp(-0.5 * (AICc - min(AICc)))),
                                                  BIC_ew = exp(-0.5 * (BIC - min(BIC))) / sum(exp(-0.5 * (BIC - min(BIC))))) 

# confidence intervals
models_list_surface_REML <- list(model_1 = surface_metrics_model_1, 
                                 model_2 = surface_metrics_model_2, 
                                 model_3 = surface_metrics_model_3, 
                                 model_4 = surface_metrics_model_4)

ci_surface_intervals <- purrr::map(models_list_surface_REML, function(x) {
  confint(x, level = 0.95, method = "Wald")
})

full_model_surface <- lme4::lmer(RST ~ Sa_scaled + S10z_scaled + Ssk + 
  Sdr_scaled + Sbi + Std_scaled + 
  Stdi + Sfd + Srwi + (1|site_1), 
  data = surface_metrics, 
  REML = FALSE, na.action = "na.fail")

model_dredge_surface <- MuMIn::dredge(full_model_surface)
head(model_dredge_surface)

#### Patch metrics ####

# import data landscape metrics
landscape_metrics <- readr::read_rds("data/Output/landscape_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add rst value
landscape_metrics <- dplyr::left_join(landscape_metrics, rst, 
                                      by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
landscape_metrics_lndscp <- dplyr::filter(landscape_metrics, 
                                          level == "landscape") %>%
  dplyr::select(site_a, site_b, metric, value, RST, euclidean_distance)

# reshape to wide format
landscape_metrics_lndscp <- tidyr::spread(landscape_metrics_lndscp, 
                                          metric, value)

# remove pr and rpr (pr = 3 and rpr = 100% for all clips)
landscape_metrics_lndscp <- dplyr::select(landscape_metrics_lndscp, 
                                          -rpr, -pr)

# Create Zl and ZZ matrix
Zl_landscape <- lapply(c("site_a","site_b"), function(x) {Matrix::fac2sparse(landscape_metrics_lndscp[[x]], "d", drop = FALSE)})

ZZ_landscape <- Reduce("+", Zl_landscape[-1], Zl_landscape[[1]])

# fitting model (not sure what the problem is...)
landscapemetrics_model <- lme4::lFormula(RST ~ ai + area_mn + cai_mn + condent + 
                                           contag + core_mn + division + ed + 
                                           ent + iji + joinent + lpi + lsi + mesh +
                                           mutinf + np + pd + pladj + prd +
                                           shdi + shei + siei + split + 
                                           ta + te + (1|site_a), 
                                        data = landscape_metrics_lndscp, REML = TRUE)

# Checking for variance inflation multicollinearity
dplyr::select(landscape_metrics_lndscp, 
              ai, area_mn, cai_mn, condent, contag, core_mn, division, ed, 
              ent, iji, joinent, lpi, lsi, mesh, mutinf, np, pd, pladj, 
              prd, shdi, shei, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)

# removing the metric with the highest value subsequently ending up with ones below 10
dplyr::select(landscape_metrics_lndscp, 
              cai_mn, core_mn, iji, mesh, pd, prd, split) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)

# fitting model (different scales)
landscapemetrics_model <- lme4::lFormula(RST ~ cai_mn + core_mn + iji + mesh + 
                                           pd + prd + split + (1|site_a), 
                                         data = landscape_metrics_lndscp, REML = TRUE)

# rescale metrics
landscape_metrics_lndscp <- dplyr::mutate(landscape_metrics_lndscp, 
                                          cai_mn_scaled = as.numeric(scale(cai_mn)), 
                                          core_mn_scaled = as.numeric(scale(core_mn)),
                                          iji_scaled = as.numeric(scale(iji)),
                                          mesh_scaled = as.numeric(scale(mesh)),
                                          pd_scaled = as.numeric(scale(pd)),
                                          prd_scaled = as.numeric(scale(prd)),
                                          split_scaled = as.numeric(scale(split)))

# update model
landscapemetrics_model_1 <- modular_function(RST ~ cai_mn_scaled + core_mn_scaled + 
                                             iji_scaled + mesh_scaled + pd_scaled + 
                                             prd_scaled +  split_scaled + (1|site_a), 
                                             data = landscape_metrics_lndscp,
                                             ZZ = ZZ_landscape, 
                                             REML = TRUE)

# look at model summary
summary(landscapemetrics_model_1)

# plot model results
ggplot() + 
  geom_point(aes(x = fitted(landscapemetrics_model_1), 
                 y = residuals(landscapemetrics_model_1)), 
             pch = 1, size = 2.5) + 
  geom_hline(yintercept = 0) + 
  labs(x = "fitted values", y = "residuals") +
  theme_bw()

# update model
landscapemetrics_model_1_nr <- modular_function(RST ~ cai_mn_scaled + core_mn_scaled + 
                                                  iji_scaled + mesh_scaled + pd_scaled + 
                                                  prd_scaled +  split_scaled + (1|site_a), 
                                                data = landscape_metrics_lndscp,
                                                ZZ = ZZ_landscape, 
                                                REML = FALSE)

landscapemetrics_model_2_nr <- modular_function(RST ~ iji_scaled + mesh_scaled + split_scaled + (1|site_a), 
                                                data = landscape_metrics_lndscp,
                                                ZZ = ZZ_landscape, 
                                                REML = FALSE)

landscapemetrics_model_3_nr <- modular_function(RST ~ cai_mn_scaled + core_mn_scaled +(1|site_a), 
                                                data = landscape_metrics_lndscp,
                                                ZZ = ZZ_landscape, 
                                                REML = FALSE)

landscapemetrics_model_4_nr <- modular_function(RST ~  pd_scaled + prd_scaled + (1|site_a), 
                                                data = landscape_metrics_lndscp,
                                                ZZ = ZZ_landscape, 
                                                REML = FALSE)

# update model
landscapemetrics_model_1 <- modular_function(RST ~ cai_mn_scaled + core_mn_scaled + 
                                               iji_scaled + mesh_scaled + pd_scaled + 
                                               prd_scaled +  split_scaled + (1|site_a), 
                                             data = landscape_metrics_lndscp,
                                             ZZ = ZZ_landscape, 
                                             REML = TRUE)

landscapemetrics_model_2 <- modular_function(RST ~ iji_scaled + mesh_scaled + split_scaled + (1|site_a), 
                                             data = landscape_metrics_lndscp,
                                             ZZ = ZZ_landscape, 
                                             REML = TRUE)

landscapemetrics_model_3 <- modular_function(RST ~ cai_mn_scaled + core_mn_scaled +(1|site_a), 
                                             data = landscape_metrics_lndscp,
                                             ZZ = ZZ_landscape, 
                                             REML = TRUE)

landscapemetrics_model_4 <- modular_function(RST ~  pd_scaled + prd_scaled + (1|site_a), 
                                             data = landscape_metrics_lndscp,
                                             ZZ = ZZ_landscape, 
                                             REML = TRUE)

# put all models in list
models_list_landscape_nr <- list(model_1 = landscapemetrics_model_1_nr, 
                                 model_2 = landscapemetrics_model_2_nr, 
                                 model_3 = landscapemetrics_model_3_nr, 
                                 model_4 = landscapemetrics_model_4_nr)

# calculating AIC and BIC
information_criterion_landscape_nr <- purrr::map_dfr(models_list_landscape_nr, function(x) {
  data.frame(AIC = AIC(x), 
             BIC = BIC(x),
             log_like = attr(logLik(x), "df"))
}, .id = "model")

# calculating corrected AICc and BICc
information_criterion_landscape_nr <- dplyr::mutate(information_criterion_landscape_nr, 
                                                    AICc = AIC + 2 * log_like * 
                                                      (log_like + 1) / (48 - log_like - 1),
                                                    AICc_min = exp(-0.5 * (AICc - min(AICc))) / sum(exp(-0.5 * (AICc - min(AICc)))),
                                                    BIC_ew = exp(-0.5 * (BIC - min(BIC))) / sum(exp(-0.5 * (BIC - min(BIC))))) 

# confidence intervals
models_list_landscape_REML <- list(model_1 = landscapemetrics_model_1, 
                                   model_2 = landscapemetrics_model_2, 
                                   model_3 = landscapemetrics_model_3, 
                                   model_4 = landscapemetrics_model_4)

ci_surface_intervals <- purrr::map(models_list_landscape_REML, function(x) {
  confint(x, level = 0.95, method = "Wald")
})


# USE MuMIn::dredge()

# fit full model
full_model <- lme4::lmer(RST ~ cai_mn_scaled + core_mn_scaled + 
  iji_scaled + mesh_scaled + pd_scaled + 
  prd_scaled +  split_scaled +
  (1|site_a), 
  data = landscape_metrics_lndscp,REML = FALSE, na.action = "na.fail")

# dredge models
model_dredge <- MuMIn::dredge(full_model)

head(model_dredge)

# only models with a delta < 2 (rule of thumb)
subset(model_dredge, delta < 2)
