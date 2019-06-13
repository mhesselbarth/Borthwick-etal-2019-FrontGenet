# model optimization function
modular_function <- function(variables, data, REML = TRUE,
                             ZZ = NULL, call = NULL) {
  
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
  
  # getting AIC and r2
  # should we use extractAIC here ?
  information_criterion <- purrr::map_dfr(model, function(x) {
    
    ic <- tibble::tibble(AIC = AIC(x),
                         BIC = BIC(x),
                         log_like = attr(logLik(x), "df"))
    
    r2 <- MuMIn::r.squaredGLMM(x)
    
    ic <- dplyr::mutate(ic, 
                        r2_marginal = r2[1, 1], 
                        r2_conditional = r2[1, 2])
  }, .id = "model")
  
  # calculating corrected AICc and BICc and weights
  information_criterion <- dplyr::mutate(information_criterion,
                                         AICc = AIC + 2 * log_like *
                                           (log_like + 1) / (n - log_like - 1),
                                         AICc_ew = exp(-0.5 * (AICc - min(AICc))) /
                                           sum(exp(-0.5 * (AICc - min(AICc)))),
                                         BIC_ew = exp(-0.5 * (BIC - min(BIC))) /
                                           sum(exp(-0.5 * (BIC - min(BIC)))))
  
  return(information_criterion)
}

create_ZZ <- function(data, col_names = c("site_1", "site_2")) {
  
  Zl <- lapply(col_names, function(x) {Matrix::fac2sparse(data[[x]], "d", 
                                                          drop = FALSE)})
  
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  
  return(ZZ)
}