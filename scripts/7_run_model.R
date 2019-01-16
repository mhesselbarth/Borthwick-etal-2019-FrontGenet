# load libraries
library(tidyverse)
library(lme4)

#### Preprocessing data ####

# import data surface metrics
surface_metrics <- readr::read_rds("data/Output/surface_metrics.rds")

# import data landscape metrics
landscape_metrics <- readr::read_rds("data/Output/landscape_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add RST to data sets
surface_metrics <- dplyr::left_join(surface_metrics, rst, by = c("site_1", "site_2"))

landscape_metrics <- dplyr::left_join(landscape_metrics, rst, 
                                      by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
landscape_metrics_lndscp <- dplyr::filter(landscape_metrics, 
                                          level == "landscape") %>%
  dplyr::select(site_a, site_b, metric, value, RST, euclidean_distance)

# reshape to wide format
landscape_metrics_lndscp <- tidyr::spread(landscape_metrics_lndscp, 
                                          metric, value)

#### Surface metrics ####

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

dplyr::select(surface_metrics, 
              Sa_scaled,S10z_scaled,Ssk,Sku,Sdr_scaled,Sbi,Std_scaled,Stdi,Sfd,Srwi) %>% 
  as.data.frame() %>%
  usdm::vif()

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

# model optimization function
modular_function <- function(variables, data, REML = TRUE) {
  
  # parse the data and formula
  model_formula <- lme4::lFormula(variables, data = data, REML = REML)
  
  # create the deviance function to be optimized
  deviance_function <- do.call(lme4::mkLmerDevfun, model_formula)
  
  # optimize the deviance function:
  optimization <- lme4::optimizeLmer(deviance_function)
  
  # package up the results
  lme4::mkMerMod(rho = environment(deviance_function), 
                 opt = optimization, 
                 reTrms = model_formula$reTrms,
                 fr = model_formula$fr)
}

# update model
surface_metrics_model_1 <- modular_function(RST ~ Sa_scaled + S10z_scaled + Ssk + 
                                              Sdr_scaled + Sbi + Std_scaled + 
                                              Stdi + Sfd + Srwi + (1|site_1), 
                                            data = surface_metrics)

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
                                            data = surface_metrics)

surface_metrics_model_3 <- modular_function(RST ~ Ssk + Sbi + Sfd + S10z_scaled + 
                                              Sdr_scaled + Std_scaled + 
                                              Sa_scaled + (1|site_1), 
                                            data = surface_metrics)

surface_metrics_model_4 <- modular_function(RST ~ Ssk + Sfd + S10z_scaled + 
                                              Sdr_scaled + Std_scaled + (1|site_1),
                                            data = surface_metrics)

# setting REML to FALSE
surface_metrics_model_1_nr <- modular_function(RST ~ Sa_scaled + S10z_scaled + Ssk + 
                                                 Sdr_scaled + Sbi + Std_scaled + 
                                                 Stdi + Sfd + Srwi + (1|site_1), 
                                               data = surface_metrics, REML = FALSE)

surface_metrics_model_2_nr <- modular_function(RST ~ Ssk + Sbi + Stdi + Sfd + Srwi + 
                                                 S10z_scaled + Sdr_scaled + Std_scaled + 
                                                 Sa_scaled + (1|site_1), 
                                               data = surface_metrics, REML = FALSE)

surface_metrics_model_3_nr <- modular_function(RST ~ Ssk + Sbi + Sfd + S10z_scaled + 
                                                 Sdr_scaled + Std_scaled + 
                                                 Sa_scaled + (1|site_1), 
                                               data = surface_metrics, REML = FALSE)

surface_metrics_model_4_nr <- modular_function(RST ~ Ssk + Sfd + S10z_scaled + 
                                                 Sdr_scaled + Std_scaled + (1|site_1),
                                               data = surface_metrics, REML = FALSE)

# put all models in list
models_list <- list(model_1 = surface_metrics_model_1_nr, 
                    model_2 = surface_metrics_model_2_nr, 
                    model_3 = surface_metrics_model_3_nr, 
                    model_4 = surface_metrics_model_4_nr)

# calculating AIC and BIC
information_criterion <- purrr::map_dfr(models_list, function(x) {
  data.frame(AIC = AIC(x), 
             BIC = BIC(x),
             log_like = attr(logLik(x), "df"))
  }, .id = "model")

# calculating corrected AICc and BICc
information_criterion <- dplyr::mutate(information_criterion, 
                                       AICc = AIC + 2 * log_like * 
                                         (log_like + 1) / (48 - log_like - 1),
                                       AICc_min = exp(-0.5 * (AICc - min(AICc))) / sum(exp(-0.5 * (AICc - min(AICc)))),
                                       BIC_ew = exp(-0.5 * (BIC - min(BIC))) / sum(exp(-0.5 * (BIC - min(BIC))))
                                       ) 

# confidence intervals
models_list_REML <- list(model_1 = surface_metrics_model_1, 
                         model_2 = surface_metrics_model_2, 
                         model_3 = surface_metrics_model_3, 
                         model_4 = surface_metrics_model_4)

ci_intervals <- purrr::map(models_list_REML, function(x) {
  confint(x, level = 0.95, method = "Wald")
})

#### Patch metrics ####
landscapemetrics_model <- lme4::lFormula(RST ~ ai + area_mn + cai_mn + condent + 
                                          contag + core_mn + division + ed + 
                                          ent + iji + joinent + lpi + lsi + mesh +
                                          mutinf + np + pd + pladj + pr + prd +
                                          rpr + shdi + shei + siei + split + 
                                          ta + te + (1|site_a), 
                                        data = landscape_metrics_lndscp, REML = TRUE)

dplyr::select(landscape_metrics_lndscp, 
              ai, area_mn, cai_mn, condent, contag, core_mn, division, ed, 
              ent, iji, joinent, lpi, lsi, mesh, mutinf, np, pd, pladj, pr, 
              prd, rpr, shdi, shei, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)


dplyr::select(landscape_metrics_lndscp, 
              cai_mn, core_mn, iji, mesh, pd,  prd, rpr, split) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)

landscapemetrics_model <- lme4::lFormula(RST ~ cai_mn + core_mn + iji + mesh + 
                                           pd + prd + rpr + split + (1|site_a), 
                                         data = landscape_metrics_lndscp, REML = TRUE)

