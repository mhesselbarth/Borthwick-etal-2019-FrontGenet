# model optimization function
modular_function <- function(response, explanatory, random,
                             data, 
                             REML = TRUE,
                             ZZ = NULL, 
                             verbose = TRUE) {
  
  model_formula <- as.formula(paste(paste(response, 
                                          paste(explanatory, 
                                                collapse = " + "), sep = " ~ "), 
                                    random, sep = " + "))
  
  
  # parse the data and formula
  model <- lme4::lFormula(model_formula, data = data, REML = REML)
  
  if (!is.null(ZZ)) {
    
    # replace ZZ matrix
    model$reTrms$Zt <- ZZ
    
    if (verbose) {  
      message("> Replaced 'model$reTrms$Zt <- ZZ'")
    }
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
  
  return(model)
}

get_model_aic <- function(model, n) {
  
  # getting AIC and r2
  # should we use extractAIC here ?
  information_criterion <- tibble::tibble(AIC = AIC(model),
                                          BIC = BIC(model),
                                          log_like = attr(logLik(model), "df"))
  
  # calculating corrected AICc and BICc and weights
  information_criterion <- dplyr::mutate(information_criterion,
                                         AICc = AIC + 2 * log_like *
                                           (log_like + 1) / (n - log_like - 1))
  
  return(information_criterion)
}

get_model_r2 <- function(model) {
  
  r2 <- tibble::as_tibble(MuMIn::r.squaredGLMM(model))
  
  return(r2)
}

create_ZZ <- function(data, col_names = c("site_1", "site_2")) {
  
  Zl <- lapply(col_names, function(x) {Matrix::fac2sparse(data[[x]], "d", 
                                                          drop = FALSE)})
  
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  
  return(ZZ)
}