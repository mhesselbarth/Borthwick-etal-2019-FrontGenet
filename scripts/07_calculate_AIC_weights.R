#### Calculate AIC weights #### 

# before this script runs, you need to run 05_run_models.R and 06_LCP_models.R to 
# have all required data in workspace

# combine AICc into one df
df_AIC <- dplyr::bind_rows(
  dplyr::bind_cols(approach = "surface", AIC = model_dredge_surface$AIC[1]),
  dplyr::bind_cols(approach = "lsm", AIC = model_dredge_lsm$AIC[1]),
  dplyr::bind_cols(approach = "lsm_NA", AIC = model_dredge_lsm_NA$AIC[1]),
  dplyr::bind_cols(approach = "ibd", AIC = model_dredge_ibd$AIC[1]),
  dplyr::bind_cols(approach = "lcp", AIC = model_dredge_lcp$AIC[1]),
  dplyr::bind_cols(approach = "res", AIC = model_dredge_res$AIC[1]))

# calculate AIC weights
df_AIC <- dplyr::mutate(df_AIC, 
                        AIC_delta = AIC - min(AIC), 
                        AIC_ew = exp(-0.5 * AIC_delta) / sum(AIC_delta)) %>%
  dplyr::arrange(AIC)
