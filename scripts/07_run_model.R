# load libraries
library(ecodist)
library(tidyverse)
library(lme4)
library(Matrix)
library(MuMIn)
library(usdm)

# source modular functions
source("scripts/00_modular_functions.R")

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
df_surface <- readr::read_rds("data/Output/surface_metrics.rds")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add RST to data sets
df_surface <- dplyr::left_join(df_surface, rst, 
                               by = c("site_1", "site_2"))

# add distance to data sets
df_surface <- dplyr::mutate(df_surface,
                            euclidean_distance = distance_matrix[cbind(site_1, 
                                                                       site_2)]) %>% 
  dplyr::arrange(site_1, site_2)

# Checking for variance inflation multicollinearity (VIF)
# removing metrics with a VIF > 10
dplyr::select(df_surface, 
              Sa, S10z, Ssk, Sku, Sdr, Sbi, Std, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# Create ZZ matrix
ZZ_surface <- create_ZZ(data = df_surface, 
                        col_names = c("site_1", "site_2"))

# rescale metrics
df_surface <- dplyr::mutate(df_surface, 
                            Sa_scaled = as.numeric(scale(Sa, 
                                                         center = TRUE, 
                                                         scale = TRUE)),
                            S10z_scaled = as.numeric(scale(S10z,
                                                           center = TRUE, 
                                                           scale = TRUE)),
                            Ssk_scaled = as.numeric(scale(Ssk,
                                                          center = TRUE, 
                                                          scale = TRUE)),
                            Sdr_scaled = as.numeric(scale(Sdr, 
                                                          center = TRUE,
                                                          scale = TRUE)), 
                            Sbi_scaled = as.numeric(scale(Sbi, 
                                                          center = TRUE, 
                                                          scale = TRUE)), 
                            Std_scaled = as.numeric(scale(Std, 
                                                          center = TRUE,
                                                          scale = TRUE)),
                            Stdi_scaled = as.numeric(scale(Stdi,
                                                           center = TRUE,
                                                           scale = TRUE)),
                            Sfd_scaled = as.numeric(scale(Sfd, 
                                                          center = TRUE, 
                                                          scale = TRUE)), 
                            Srwi_scaled = as.numeric(scale(Srwi,
                                                           center = TRUE, 
                                                           scale = TRUE)),
                            dist_scaled = as.numeric(scale(euclidean_distance,
                                                           center = TRUE, 
                                                           scale = TRUE)))

# specify full models 
model_full_no_REML_surface <- lme4::lmer(formula = RST ~ Sa_scaled + 
                                           S10z_scaled + Ssk_scaled +
                                           Sdr_scaled + Sbi_scaled + 
                                           Std_scaled + Stdi_scaled + 
                                           Sfd_scaled + Srwi_scaled +
                                           (1|site_1), 
                                         data = df_surface, 
                                         REML = FALSE, 
                                         na.action = "na.fail")

# dredge lme model
model_dredge_surface <- MuMIn::dredge(model_full_no_REML_surface)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_surface <- MuMIn::get.models(model_dredge_surface, 
                                              subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_surface <- purrr::map_dfr(seq_along(all_combinations_surface), function(x) {
  
  message("\r> Progress: ", x, "/", length(all_combinations_surface), 
          appendLF = FALSE)
  
  model_fit <- modular_function(response = "RST", 
                   explanatory = all_combinations_surface[[x]],
                   random =  "(1|site_1)", 
                   data = df_surface, 
                   REML = FALSE, 
                   ZZ = ZZ_surface,
                   verbose = FALSE)
  
  get_model_aic(model_fit, n = 136)
}, .id = "model")

# order by AICc
model_dredge_surface <- dplyr::arrange(model_dredge_surface, 
                                       AICc) %>%
  dplyr::mutate(model = as.numeric(model))

# best models
best_model_surface <- all_combinations_surface[model_dredge_surface$model[1:3]]

# rerun best models using REML = TRUE and extract R2
r2_surface <- purrr::map_dfr(model_dredge_surface$model[1:3], function(x) {
  
  model_fit <- modular_function(response = "RST", 
                                explanatory = all_combinations_surface[[x]],
                                random =  "(1|site_1)", 
                                data = df_surface, 
                                REML = TRUE, 
                                ZZ = ZZ_surface,
                                verbose = FALSE)
  
  dplyr::bind_cols(model = x, get_model_r2(model_fit))
})

#### Patch metrics ####

# import data landscape metrics
df_lsm <- readr::read_rds("data/Output/landscape_metrics.rds")

# only landscape level metrics
df_lsm <- dplyr::filter(df_lsm, level == "landscape")

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add rst value
df_lsm <- dplyr::left_join(df_lsm, rst,
                           by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
df_lsm <- dplyr::select(df_lsm, 
                        site_a, site_b, 
                        metric, value, 
                        RST, euclidean_distance)

# reshape to wide format
df_lsm <- tidyr::spread(df_lsm, 
                        metric, value)

# remove pr and rpr (pr = 3 and rpr = 100% for all clips)
df_lsm <- dplyr::select(df_lsm, 
                        -rpr, -pr)

# Checking for variance inflation multicollinearity
# Removing the metric with the highest value subsequently ending up with ones below 10
dplyr::select(df_lsm, 
              ai, area_mn, cai_mn, condent, contag, core_mn, division, ed, 
              ent, iji, joinent, lpi, lsi, mesh, mutinf, np, pd, pladj, 
              prd, shdi, shei, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# Create Zl and ZZ matrix
ZZ_landscape <- create_ZZ(data = df_lsm, col_names = c("site_a", "site_b"))

# rescale metrics
df_lsm <- dplyr::mutate(df_lsm, 
                        pd_scaled = as.numeric(scale(pd)),
                        iji_scaled = as.numeric(scale(iji)),
                        prd_scaled = as.numeric(scale(prd)),
                        core_mn_scaled = as.numeric(scale(core_mn)),
                        cai_mn_scaled = as.numeric(scale(cai_mn)), 
                        split_scaled = as.numeric(scale(split)),
                        mesh_scaled = as.numeric(scale(mesh)),
                        dist_scaled = as.numeric(scale(euclidean_distance)))

# specify full models 
model_full_no_REML_lsm <- lme4::lmer(formula = RST ~ pd_scaled + 
                                       iji_scaled + prd_scaled + 
                                       core_mn_scaled + cai_mn_scaled + 
                                       split_scaled + mesh_scaled + 
                                       (1|site_a), 
                                     data = df_lsm, 
                                     REML = FALSE, 
                                     na.action = "na.fail")

# dredge lme model
model_dredge_lsm <- MuMIn::dredge(model_full_no_REML_lsm)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_lsm <- MuMIn::get.models(model_dredge_lsm,
                                          subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_lsm <- purrr::map_dfr(seq_along(all_combinations_lsm), function(x) {
  
  message("\r> Progress: ", x, "/", length(all_combinations_lsm), 
          appendLF = FALSE)
  
  model_fit <- modular_function(response = "RST", 
                                explanatory = all_combinations_lsm[[x]],
                                random =  "(1|site_a)", 
                                data = df_lsm, 
                                REML = FALSE, 
                                ZZ = ZZ_landscape,
                                verbose = FALSE)
  
  get_model_aic(model_fit, n = 136)
}, .id = "model")

# order by AICc
model_dredge_lsm <- dplyr::arrange(model_dredge_lsm, 
                                   AICc) %>%
  dplyr::mutate(model = as.numeric(model))

# best models
best_model_lsm <- all_combinations_lsm[model_dredge_lsm$model[1:3]]

# rerun best models using REML = TRUE and extract R2
r2_lsm <- purrr::map_dfr(model_dredge_lsm$model[1:3], function(x) {
  
  model_fit <- modular_function(response = "RST", 
                                explanatory = all_combinations_lsm[[x]],
                                random =  "(1|site_a)", 
                                data = df_lsm, 
                                REML = TRUE, 
                                ZZ = ZZ_landscape,
                                verbose = FALSE)
  
  dplyr::bind_cols(model = x, get_model_r2(model_fit))
})

#### Isolation by distance

# import RST value
rst <- readr::read_rds("data/rst.rds")

# add distance to rst
df_ibd <- dplyr::mutate(rst,
                        euclidean_distance = distance_matrix[cbind(site_1,
                                                                   site_2)],
                     dist_scaled = as.numeric(scale(euclidean_distance))) %>% 
  dplyr::arrange(site_1, site_2) 

# create ZZ matrix
Zl_dist <- lapply(c("site_1","site_2"), function(x) {
  Matrix::fac2sparse(df_ibd[[x]], "d", drop = FALSE)})

ZZ_dist <- Reduce("+", Zl_dist[-1], Zl_dist[[1]])

# fit model using ZZ matrix
model_dredge_ibd <- modular_function(response = "RST", 
                                   explanatory = "dist_scaled", 
                                   random = "(1|site_1)", 
                                   data = df_ibd,
                                   ZZ = ZZ_dist,
                                   REML = FALSE) %>%
  get_model_aic(n = 136) %>%
  dplyr::bind_cols(model = 1, .)

r2_ibd <- modular_function(response = "RST", 
                           explanatory = "dist_scaled", 
                           random = "(1|site_1)", 
                           data = df_ibd,
                           ZZ = ZZ_dist,
                           REML = TRUE) %>% 
  get_model_r2() %>%
  dplyr::bind_cols(model = 1, .)

#### Compare all three models ####
best_model_surface
best_model_lsm

model_dredge_surface[1:3, ]
model_dredge_lsm[1:3, ]
model_dredge_ibd

r2_surface
r2_lsm
r2_ibd
