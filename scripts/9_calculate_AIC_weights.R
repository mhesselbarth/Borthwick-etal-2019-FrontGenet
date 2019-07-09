#### Calculate AIC weights #### 

# before this script runs, you need to run 7_run_models.R and 8_LCP_models.R to 
# have all required data in workspace

# combine AICc into one df
df_AICc <- dplyr::bind_rows(
  dplyr::bind_cols(approach = "surface", AICc = model_dredge_surface$AICc[1]),
  dplyr::bind_cols(approach = "lsm", AICc = model_dredge_lsm$AICc[1]),
  dplyr::bind_cols(approach = "ibd", AICc = model_dredge_ibd$AICc[1]),
  dplyr::bind_cols(approach = "lcp", AICc = model_dredge_lcp$AICc[1]),
  dplyr::bind_cols(approach = "res", AICc = model_dredge_res$AICc[1]))

# calculate AIC weights
df_AICc <- dplyr::mutate(df_AICc, 
                         AICc_delta = AICc - min(AICc), 
                         AIC_ew = exp(-0.5 * AICc_delta) / sum(AICc_delta)) %>%
  dplyr::arrange(AICc)
