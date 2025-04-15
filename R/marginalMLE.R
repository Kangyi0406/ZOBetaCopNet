#' Marginal Maximum Likelihood Estimation
#'
#' Maximum likelihood estimation of marginal distribution function parameters.
#'
#' @param formula.sigma an object of class formula: a symbolic description of the model to be fitted for the beta dispersion parameter.
#' @param formula.mu an object of class formula: a symbolic description of the model to be fitted for the beta mean parameter.
#' @param formula.nu an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param formula.tau an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param df data frame of outcomes and covariates.
#'
#' @return \code{marginalMLE} returns a list of the following elements:
#' \item{coefficients}{vector of regression coefficients from the models for nu, tau, sigma, and mu.}
#' \item{fitted}{data frame of fitted values for nu, tau, sigma, and mu.}
#'
#' @export
#' 
marginalMLE = function(formula.mu = ~ 1, formula.sigma = ~ 1, formula.nu = ~ 1, formula.tau = ~ 1, df = NULL){
  if(!("data.frame" %in% class(df))) {
    stop("ERROR: df must be a data frame of covariates for the marginal regression models.")
  }
  #formula.mu = y~1; formula.sigma = ~ 1; formula.p = ~ 1; formula.tau = ~ 1
  # make inputted formula of class "formula"
  formula.mu = stats::as.formula(formula.mu)
  formula.sigma = stats::as.formula(formula.sigma)
  formula.nu = stats::as.formula(formula.nu)
  formula.tau = stats::as.formula(formula.tau)
  
  # fit ZI-Beta regression model
  # model = suppressWarnings( # suppress function warnings
  #   gamlss::gamlss(formula = formula.mu, sigma.formula = formula.sigma, nu.formula = formula.p, formula.tau = formula.tau,
  #                  family = gamlss.dist::BEINF, data = df,
  #                  control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  #   # gamlss.inf::gamlssInf0to1(Y = df[,"y"], mu.formula = formula.mu, sigma.formula = formula.sigma, xi0.formula = formula.p, xi1.formula = formula.tau,
  #   #                family = gamlss.dist::BE, data = df,
  #   #                control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  #   
  # )
  model = gamlss::gamlss(formula = formula.mu, sigma.formula = formula.sigma, nu.formula = formula.nu, formula.tau = formula.tau,
                         family = gamlss.dist::BEINF, data = df,
                         control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  # gamlss.inf::gamlssInf0to1(Y = df[,"y"], mu.formula = formula.mu, sigma.formula = formula.sigma, xi0.formula = formula.nu, xi1.formula = formula.tau,
  #                family = gamlss.dist::BE, data = df,
  #                control = gamlss::gamlss.control(n.cyc = 40, trace=FALSE))
  
  # Capture warnings
  warn <- warnings()
  
  # Count the number of warnings
  num_warnings <- length(warn)
  # check convergence
  convergenceFail = model$converged==FALSE
  
  # if model failed to converge refit
  if(convergenceFail){
    
    message("Model failed to converge.")
    
    
    model = gamlss::refit(model)
    
    if(model$converged==FALSE){
      message("Model refit and failed to converged.")
    } else {
      message("Model successfully refit.")
    }
  }
  
  
  
  # combine estimated regression coefficients
  zibrCoef = unlist(gamlss::coefAll(model))
  names(zibrCoef) = janitor::make_clean_names(names(zibrCoef))
  
  # estimated zero probability, mean, and dispersion
  mu = stats::fitted(model,"mu"); sigma = stats::fitted(model,"sigma"); nu = stats::fitted(model,"nu"); tau = stats::fitted(model,"tau")
  alpha = mu*(1-sigma^2)/sigma^2
  beta = (1-mu)*(1-sigma^2)/sigma^2
  p = nu/(1+nu+tau)
  q = tau/(1+nu+tau)
  
  
  
  # output
  out = list("coefficients" = zibrCoef,
             "fitted" = data.frame(q=q, p=p, alpha=alpha, beta=beta),
             "num_warnings" = num_warnings
  )
  
  return(out)
}

