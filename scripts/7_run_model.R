# load libraries
library(ecodist)
library(tidyverse)
library(lme4)
library(Matrix)
library(MuMIn)
library(usdm)

# model optimization function
modular_function <- function(variables, data, REML = TRUE, ZZ = NULL, call = NULL) {

  # parse the data and formula
  model <- lme4::lFormula(variables, data = data, REML = REML)
  
  if (!is.null(ZZ)) {
    
    # replace ZZ matrix
    model$reTrms$Zt <- ZZ
    
    message("> Replaced 'model$reTrms$Zt <- ZZ'")
  }

  # create the deviance function to be optimized
  deviance_function <- do.call(lme4::mkLmerDevfun, model)

  # optimize the deviance function:
  optimization <- lme4::optimizeLmer(deviance_function)

  # package up the results
  model <- lme4::mkMerMod(rho = environment(deviance_function),
                          opt = optimization,
                          reTrms = model$reTrms,
                          fr = model$fr)
  
  if (!is.null(call)) {
    
    # replace call for dredging
    model@call <- call
    
    message("> Replaced 'model@call <- call'")
    
  }
  
  return(model)
}

get_model_info <- function(model, n) {
  
  information_criterion <- tibble::tibble(AIC = AIC(model),
                                          BIC = BIC(model),
                                          log_like = attr(logLik(model), "df"))

  # calculating corrected AICc and BICc
  # not sure if weigths are correct at the moment
  information_criterion <- dplyr::mutate(information_criterion,
                                         AICc = AIC + 2 * log_like *
                                           (log_like + 1) / (n - log_like - 1),
                                         AICc_ew = exp(-0.5 * (AICc - min(AICc))) /
                                           sum(exp(-0.5 * (AICc - min(AICc)))),
                                         BIC_ew = exp(-0.5 * (BIC - min(BIC))) /
                                           sum(exp(-0.5 * (BIC - min(BIC)))))
  
  r2 <- MuMIn::r.squaredGLMM(model)
  
  information_criterion <- dplyr::mutate(information_criterion, 
                                         r2_marginal = r2[1, 1], 
                                         r2_conditional = r2[1, 2])

  return(information_criterion)
}

#### Import distance

# Import sample points
sample_points <- getwd() %>%
  paste0("/data/GIS/SSR_17_sites.shp") %>%
  raster::shapefile() %>%
  tibble::as_tibble()

# Calculate distance between each point
distance_matrix <- sample_points %>%
  dplyr::select(E, N) %>%
  dist(diag = TRUE, upper = TRUE) %>%
  as.matrix()

#### Surface metrics ####

# import data surface metrics
surface_metrics <- readr::read_rds("data/Output/surface_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add RST to data sets
surface_metrics <- dplyr::left_join(surface_metrics, rst, 
                                    by = c("site_1", "site_2"))

# add distance to data sets
surface_metrics <- dplyr::mutate(surface_metrics,
                                 euclidean_distance = distance_matrix[cbind(site_1, 
                                                                            site_2)]) %>% 
  dplyr::arrange(site_1, site_2)

# Checking for variance inflation multicollinearity (VIF)
dplyr::select(surface_metrics, 
              Sa, S10z, Ssk, Sku, Sdr, Sbi, Std, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vif()

# look at VIF, removed Sku
dplyr::select(surface_metrics, 
              Sa, S10z, Ssk, Sdr, Sbi, Std, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vif()

# Create Zl and ZZ matrix
Zl_surface <- lapply(c("site_1","site_2"), function(x) {Matrix::fac2sparse(surface_metrics[[x]], "d", drop = FALSE)})

ZZ_surface <- Reduce("+", Zl_surface[-1], Zl_surface[[1]])

# rescale metrics
surface_metrics <- dplyr::mutate(surface_metrics, 
                                 Sa_scaled = as.numeric(scale(Sa, center = TRUE, scale = TRUE)),
                                 S10z_scaled = as.numeric(scale(S10z, center = TRUE, scale = TRUE)),
                                 Ssk_scaled = as.numeric(scale(Ssk, center = TRUE, scale = TRUE)),
                                 Sdr_scaled = as.numeric(scale(Sdr, center = TRUE, scale = TRUE)), 
                                 Sbi_scaled = as.numeric(scale(Sbi, center = TRUE, scale = TRUE)), 
                                 Std_scaled = as.numeric(scale(Std, center = TRUE, scale = TRUE)),
                                 Stdi_scaled = as.numeric(scale(Stdi, center = TRUE, scale = TRUE)),
                                 Sfd_scaled = as.numeric(scale(Sfd, center = TRUE, scale = TRUE)), 
                                 Srwi_scaled = as.numeric(scale(Srwi, center = TRUE, scale = TRUE)),
                                 dist_scaled = as.numeric(scale(euclidean_distance, center = TRUE, scale = TRUE)))

# specify full models 
surface_metrics_model_full_no_REML <- lme4::lmer(formula = RST ~ Sa_scaled + S10z_scaled + Ssk_scaled + 
                                                   Sdr_scaled + Sbi_scaled + Std_scaled + Stdi_scaled + 
                                                   Sfd_scaled + Srwi_scaled + dist_scaled + (1|site_1), 
                                                 data = surface_metrics, 
                                                 REML = FALSE, na.action = "na.fail")

# create call object (needed for model dredging)
fun_call_no_REML <- surface_metrics_model_full_no_REML@call

# run model
surface_metrics_model_full_no_REML <- modular_function(variables = RST ~ Sa_scaled + S10z_scaled + 
                                                         Ssk_scaled + Sdr_scaled + Sbi_scaled + 
                                                         Std_scaled + Stdi_scaled + Sfd_scaled + 
                                                         Srwi_scaled + dist_scaled + (1|site_1), 
                                                       data = surface_metrics, 
                                                       REML = FALSE, 
                                                       ZZ = ZZ_surface,
                                                       call = fun_call_no_REML)

# # look at model summary
# summary(surface_metrics_model_full)

# # plot model results
# ggplot() + 
#   geom_point(aes(x = fitted(surface_metrics_model_full), 
#                  y = residuals(surface_metrics_model_full)), 
#              pch = 1, size = 2.5) + 
#   geom_hline(yintercept = 0) + 
#   labs(x = "fitted values", y = "residuals") +
#   theme_bw()

# dredge model - don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_surface_no_REML <- MuMIn::dredge(surface_metrics_model_full_no_REML) %>% 
  head(n = 5)

# fit best model
surface_metrics_best_REML <- modular_function(variables = RST ~ dist_scaled + S10z_scaled + 
                                                Sfd_scaled + Stdi_scaled + (1|site_1),
                                              data = surface_metrics,
                                              ZZ = ZZ_surface,
                                              REML = TRUE)

# confidence intervals
ci_surface_intervals <- confint(surface_metrics_best_REML, level = 0.95, 
                                method = "Wald")

info_surface_metrics_best_REML <- get_model_info(model = surface_metrics_best_REML, 
                                                 n = 136)

#### Patch metrics ####

# import data landscape metrics
landscape_metrics <- readr::read_rds("data/Output/landscape_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add rst value
landscape_metrics <- dplyr::left_join(landscape_metrics, rst, 
                                      by = c("site_a" = "site_1", "site_b" = "site_2"))

# add distance
landscape_metrics <- dplyr::mutate(landscape_metrics,
                                   euclidean_distance = distance_matrix[cbind(site_a, site_b)]) %>% 
  dplyr::arrange(site_a, site_b)

# only metrics on landscape level and needed cols
landscape_metrics <- dplyr::filter(landscape_metrics, 
                                   level == "landscape") %>%
  dplyr::select(site_a, site_b, metric, value, RST, euclidean_distance)

# reshape to wide format
landscape_metrics <- tidyr::spread(landscape_metrics, 
                                   metric, value)

# remove pr and rpr (pr = 3 and rpr = 100% for all clips)
landscape_metrics <- dplyr::select(landscape_metrics, 
                                   -rpr, -pr)

# Checking for variance inflation multicollinearity
dplyr::select(landscape_metrics, 
              ai, area_mn, cai_mn, condent, contag, core_mn, division, ed, 
              ent, iji, joinent, lpi, lsi, mesh, mutinf, np, pd, pladj, 
              prd, shdi, shei, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)

# removing the metric with the highest value subsequently ending up with ones below 10
dplyr::select(landscape_metrics, 
              cai_mn, core_mn, iji, mesh, pd, prd, split) %>% 
  as.data.frame() %>%
  usdm::vif() %>% 
  dplyr::arrange(-VIF)

# Create Zl and ZZ matrix
Zl_landscape <- lapply(c("site_a","site_b"), function(x) {Matrix::fac2sparse(landscape_metrics[[x]], "d", drop = FALSE)})

ZZ_landscape <- Reduce("+", Zl_landscape[-1], Zl_landscape[[1]])

# rescale metrics
landscape_metrics <- dplyr::mutate(landscape_metrics, 
                                   pd_scaled = as.numeric(scale(pd)),
                                   iji_scaled = as.numeric(scale(iji)),
                                   prd_scaled = as.numeric(scale(prd)),
                                   core_mn_scaled = as.numeric(scale(core_mn)),
                                   cai_mn_scaled = as.numeric(scale(cai_mn)), 
                                   split_scaled = as.numeric(scale(split)),
                                   mesh_scaled = as.numeric(scale(mesh)),
                                   dist_scaled = as.numeric(scale(euclidean_distance)))


# specify full models 
landscape_metrics_model_full_no_REML <- lme4::lmer(formula = RST ~ pd_scaled + iji_scaled +
                                                     prd_scaled + core_mn_scaled + 
                                                     cai_mn_scaled + split_scaled +
                                                     mesh_scaled + dist_scaled + 
                                                     (1|site_a), 
                                                   data = landscape_metrics, 
                                                   REML = FALSE, 
                                                   na.action = "na.fail")

# create call object (needed for model dredging)
fun_call_no_REML <- landscape_metrics_model_full_no_REML@call

# run model
landscape_metrics_model_full_no_REML <- modular_function(variables = RST ~ pd_scaled + iji_scaled + 
                                                           prd_scaled + core_mn_scaled + 
                                                           cai_mn_scaled + split_scaled +
                                                           mesh_scaled + dist_scaled + 
                                                           (1|site_a), 
                                                         data = landscape_metrics,
                                                         REML = FALSE,
                                                         ZZ = ZZ_landscape, 
                                                         call = fun_call_no_REML)

# # look at model summary
# summary(landscape_metrics_model_full)

# # plot model results
# ggplot() + 
#   geom_point(aes(x = fitted(landscape_metrics_model_full), 
#                  y = residuals(landscape_metrics_model_full)), 
#              pch = 1, size = 2.5) + 
#   geom_hline(yintercept = 0) + 
#   labs(x = "fitted values", y = "residuals") +
#   theme_bw()

# dredge model - don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_landscape_no_REML <- MuMIn::dredge(landscape_metrics_model_full_no_REML) %>%
  head(n = 5)

# fit best model
landscape_metrics_best_REML <- modular_function(variables = RST ~ dist_scaled +
                                                  mesh_scaled + (1|site_a),
                                                data = landscape_metrics,
                                                ZZ = ZZ_landscape, 
                                                REML = TRUE)

# # confidence intervals
ci_landscape_intervals <- confint(landscape_metrics_best_REML, 
                                  level = 0.95, method = "Wald")

info_landscape_metrics_best_REML <- get_model_info(model = landscape_metrics_best_REML, 
                                                   n = 136)

#### Isolation by distance

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add distance to rst
rst <- dplyr::mutate(rst, 
                     euclidean_distance = distance_matrix[cbind(site_1, 
                                                                site_2)],
                     dist_scaled = as.numeric(scale(euclidean_distance))) %>% 
  dplyr::arrange(site_1, site_2) 

# create ZZ matrix
Zl_dist <- lapply(c("site_1","site_2"), function(x) {Matrix::fac2sparse(rst[[x]], "d", drop = FALSE)})

ZZ_dist <- Reduce("+", Zl_dist[-1], Zl_dist[[1]])

# fit model using ZZ matrix
ibd_model_REML <- modular_function(variables = RST ~ dist_scaled + (1|site_1), 
                                   data = rst,
                                   ZZ = ZZ_dist,
                                   REML = TRUE)

# get model information 
get_model_info(ibd_model_REML, n = 136)
