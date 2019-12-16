# Re-doing models with Proportion of Shared Alleles (Dps) # 

# load libraries
library(ecodist)
library(tidyverse)
library(lme4)
library(Matrix)
library(MuMIn)
library(usdm)

# source modular functions
source("scripts/00_modular_functions.R")

#### Import distance and RST

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

# import Dps value
dps <- readr::read_csv("data/Prop_shared_allele.csv")

# separating elnads column first
dps <- tidyr::separate(data = dps, col = "Landscape", sep = "_",
                       into = c("Elands", "site_1", "site_2")) %>% 
  dplyr::select(-Elands) %>% 
  dplyr::mutate(site_1 = as.numeric(site_1),
                site_2 = as.numeric(site_2)) %>% 
  purrr::set_names(c("site_1", "site_2", "DPS"))

#### Surface metrics ####

# import data surface metrics
df_surface <- readr::read_rds("data/Output/surface_metrics.rds")

# add distance to data sets
df_surface <- dplyr::mutate(df_surface,
                            euclidean_distance = distance_matrix[cbind(site_1, 
                                                                       site_2)])

# add DPS to data sets
df_surface <- dplyr::left_join(df_surface, dps, by = c("site_1", "site_2"))

# Checking for variance inflation multicollinearity (VIF)
# removing metrics with a VIF > 10
dplyr::select(df_surface, 
              Sa, S10z, Ssk, Sku, Sdr, Sbi, Std, Stdi, Sfd, Srwi) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1        Sa 3.612458
# 2      S10z 6.295996
# 3       Ssk 2.810541
# 4       Sdr 1.282746
# 5       Sbi 2.240350
# 6       Std 3.668604
# 7      Stdi 4.620818
# 8       Sfd 4.823973
# 9      Srwi 6.017562

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
model_full_no_REML_surface <- lme4::lmer(formula = DPS ~ Sa_scaled + 
                                           S10z_scaled + Ssk_scaled +
                                           Sdr_scaled + Sbi_scaled + 
                                           Std_scaled + Stdi_scaled + 
                                           Sfd_scaled + Srwi_scaled +
                                           (1|site_1), 
                                         data = df_surface, 
                                         REML = FALSE, 
                                         na.action = "na.fail")

# dredge lme model
model_dredge_surface <- MuMIn::dredge(model_full_no_REML_surface, trace = 2)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_surface <- MuMIn::get.models(model_dredge_surface, 
                                              subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_surface <- purrr::map_dfr(seq_along(all_combinations_surface), function(x) {
  
  message("\r> Progress: ", x, "/", length(all_combinations_surface), 
          appendLF = FALSE)
  
  model_fit <- modular_function(response = "DPS", 
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
  
  model_fit <- modular_function(response = "DPS", 
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

# add rst value
df_lsm <- dplyr::left_join(df_lsm, dps,
                           by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
df_lsm <- dplyr::select(df_lsm, 
                        site_a, site_b, 
                        metric, value, 
                        DPS, euclidean_distance)

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
              prd, shdi, shei, sidi, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1    cai_mn 4.281565
# 2   core_mn 5.631005
# 3       iji 6.681443
# 4      mesh 1.890937
# 5        pd 7.144912
# 6       prd 5.978124
# 7     split 2.235872

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
model_full_no_REML_lsm <- lme4::lmer(formula = DPS ~ pd_scaled + 
                                       iji_scaled + prd_scaled + 
                                       core_mn_scaled + cai_mn_scaled + 
                                       split_scaled + mesh_scaled + 
                                       (1|site_a), 
                                     data = df_lsm, 
                                     REML = FALSE, 
                                     na.action = "na.fail")

# dredge lme model
model_dredge_lsm <- MuMIn::dredge(model_full_no_REML_lsm, trace = 2)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_lsm <- MuMIn::get.models(model_dredge_lsm,
                                          subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_lsm <- purrr::map_dfr(seq_along(all_combinations_lsm), function(x) {
  
  message("\r> Progress: ", x, "/", length(all_combinations_lsm), 
          appendLF = FALSE)
  
  model_fit <- modular_function(response = "DPS", 
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
  
  model_fit <- modular_function(response = "DPS", 
                                explanatory = all_combinations_lsm[[x]],
                                random =  "(1|site_a)", 
                                data = df_lsm, 
                                REML = TRUE, 
                                ZZ = ZZ_landscape,
                                verbose = FALSE)
  
  dplyr::bind_cols(model = x, get_model_r2(model_fit))
})

#### Patch metrics ####

# import data landscape metrics
df_lsm <- readr::read_rds("data/Output/landscape_metrics.rds")

# only landscape level metrics
df_lsm <- dplyr::filter(df_lsm, level == "landscape")

# add rst value
df_lsm <- dplyr::left_join(df_lsm, dps,
                           by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
df_lsm <- dplyr::select(df_lsm, 
                        site_a, site_b, 
                        metric, value, 
                        DPS, euclidean_distance)

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
              prd, shdi, shei, sidi, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1    cai_mn 4.281565
# 2   core_mn 5.631005
# 3       iji 6.681443
# 4      mesh 1.890937
# 5        pd 7.144912
# 6       prd 5.978124
# 7     split 2.235872

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
model_full_no_REML_lsm <- lme4::lmer(formula = DPS ~ pd_scaled + 
                                       iji_scaled + prd_scaled + 
                                       core_mn_scaled + cai_mn_scaled + 
                                       split_scaled + mesh_scaled + 
                                       (1|site_a), 
                                     data = df_lsm, 
                                     REML = FALSE, 
                                     na.action = "na.fail")

# dredge lme model
model_dredge_lsm <- MuMIn::dredge(model_full_no_REML_lsm, trace = 2)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_lsm <- MuMIn::get.models(model_dredge_lsm,
                                          subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_lsm <- purrr::map_dfr(seq_along(all_combinations_lsm), function(x) {
  
  message("\r> Progress: ", x, "/", length(all_combinations_lsm), 
          appendLF = FALSE)
  
  model_fit <- modular_function(response = "DPS", 
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
  
  model_fit <- modular_function(response = "DPS", 
                                explanatory = all_combinations_lsm[[x]],
                                random =  "(1|site_a)", 
                                data = df_lsm, 
                                REML = TRUE, 
                                ZZ = ZZ_landscape,
                                verbose = FALSE)
  
  dplyr::bind_cols(model = x, get_model_r2(model_fit))
})

#### Patch metrics NA ####

# import data landscape metrics
df_lsm_NA <- readr::read_rds("data/Output/landscape_metrics_NA.rds")

# only landscape level metrics
df_lsm_NA <- dplyr::filter(df_lsm_NA, level == "landscape")

# add rst value
df_lsm_NA <- dplyr::left_join(df_lsm_NA, dps,
                              by = c("site_a" = "site_1", "site_b" = "site_2"))

# only metrics on landscape level and needed cols
df_lsm_NA <- dplyr::select(df_lsm_NA, 
                           site_a, site_b, 
                           metric, value, 
                           DPS, euclidean_distance)

# reshape to wide format
df_lsm_NA <- tidyr::spread(df_lsm_NA, 
                           metric, value)

# remove pr and rpr (pr = 3 and rpr = 100% for all clips) and IJI (NA for < 3 classes)
df_lsm_NA <- dplyr::select(df_lsm_NA, 
                           -rpr, -pr, -iji)

# Checking for variance inflation multicollinearity
# Removing the metric with the highest value subsequently ending up with ones below 10
dplyr::select(df_lsm_NA, 
              ai, area_mn, cai_mn, condent, contag, core_mn, division, ed, 
              ent, joinent, lpi, lsi, mesh, mutinf, np, pd, pladj, 
              prd, shdi, shei, sidi, siei, split, ta, te) %>% 
  as.data.frame() %>%
  usdm::vifstep(th = 10)

# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1   area_mn 5.698623
# 2    cai_mn 2.016140
# 3        ed 5.345607
# 4      mesh 5.059306
# 5        np 5.821331
# 6        pd 6.597325
# 7       prd 3.000816
# 8     split 2.578870

# Create Zl and ZZ matrix
ZZ_landscape_NA <- create_ZZ(data = df_lsm_NA, col_names = c("site_a", "site_b"))

# rescale metrics
df_lsm_NA <- dplyr::mutate(df_lsm_NA, 
                           area_mn_scaled = as.numeric(scale(area_mn)),
                           cai_mn_scaled = as.numeric(scale(cai_mn)), 
                           ed_scaled = as.numeric(scale(ed)),
                           mesh_scaled = as.numeric(scale(mesh)), 
                           np_scaled = as.numeric(scale(np)),
                           pd_scaled = as.numeric(scale(pd)),
                           prd_scaled = as.numeric(scale(prd)),
                           split_scaled = as.numeric(scale(split)),
                           dist_scaled = as.numeric(scale(euclidean_distance)))

# specify full models 
model_full_no_REML_lsm_NA <- lme4::lmer(formula = DPS ~ area_mn_scaled +
                                          cai_mn_scaled + ed_scaled + mesh_scaled +
                                          np_scaled + pd_scaled +  prd_scaled + 
                                          split_scaled + (1|site_a), 
                                        data = df_lsm_NA, 
                                        REML = FALSE, 
                                        na.action = "na.fail")

# dredge lme model
model_dredge_lsm_NA <- MuMIn::dredge(model_full_no_REML_lsm_NA, trace = 2)

# get al function calls (i.e. explanatory variable combinations) of dredge
all_combinations_lsm_NA <- MuMIn::get.models(model_dredge_lsm_NA,
                                             subset = NA) %>%
  purrr::map(function(x) names(x@frame)[-c(1, length(names(x@frame)))])

# fit models using modular_functions and all combinations
# don't use REML here: https://tinyurl.com/yagpx4bv
model_dredge_lsm_NA <- purrr::map_dfr(seq_along(all_combinations_lsm_NA), 
                                      function(x) {
                                        
                                        message("\r> Progress: ", x, "/", length(all_combinations_lsm_NA), 
                                                appendLF = FALSE)
                                        
                                        model_fit <- modular_function(response = "DPS", 
                                                                      explanatory = all_combinations_lsm_NA[[x]],
                                                                      random =  "(1|site_a)", 
                                                                      data = df_lsm_NA, 
                                                                      REML = FALSE, 
                                                                      ZZ = ZZ_landscape_NA,
                                                                      verbose = FALSE)
                                        
                                        get_model_aic(model_fit, n = 136)
                                      }, .id = "model")

# order by AICc
model_dredge_lsm_NA <- dplyr::arrange(model_dredge_lsm_NA,
                                      AICc) %>%
  dplyr::mutate(model = as.numeric(model))

# best models
best_model_lsm_NA <- all_combinations_lsm_NA[model_dredge_lsm_NA$model[1:3]]

# rerun best models using REML = TRUE and extract R2
r2_lsm_NA <- purrr::map_dfr(model_dredge_lsm_NA$model[1:3], function(x) {
  
  model_fit <- modular_function(response = "DPS", 
                                explanatory = all_combinations_lsm_NA[[x]],
                                random =  "(1|site_a)", 
                                data = df_lsm_NA, 
                                REML = TRUE, 
                                ZZ = ZZ_landscape_NA,
                                verbose = FALSE)
  
  dplyr::bind_cols(model = x, get_model_r2(model_fit))
})

#### Isolation by distance
# add distance to dps
df_ibd <- dplyr::mutate(dps,
                        euclidean_distance = distance_matrix[cbind(site_1,
                                                                   site_2)],
                        dist_scaled = as.numeric(scale(euclidean_distance))) %>% 
  dplyr::arrange(site_1, site_2) 

# create ZZ matrix
Zl_dist <- lapply(c("site_1","site_2"), function(x) {
  Matrix::fac2sparse(df_ibd[[x]], "d", drop = FALSE)})

ZZ_dist <- Reduce("+", Zl_dist[-1], Zl_dist[[1]])

# fit model using ZZ matrix
model_dredge_ibd <- modular_function(response = "DPS", 
                                     explanatory = "dist_scaled", 
                                     random = "(1|site_1)", 
                                     data = df_ibd,
                                     ZZ = ZZ_dist,
                                     REML = FALSE) %>%
  get_model_aic(n = 136) %>%
  dplyr::bind_cols(model = 1, .)

r2_ibd <- modular_function(response = "DPS", 
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
best_model_lsm_NA

model_dredge_surface[1:3, ]
model_dredge_lsm[1:3, ]
model_dredge_lsm_NA[1:3, ]
model_dredge_ibd

r2_surface
r2_lsm
r2_lsm_NA
r2_ibd
