#' Simulate the data pairs follows zero-one inflation data with frank copula association

#' @param n size of x data frame with two columns for bivariate pair (x1, x2).
#' @param theta copula dependence parameter.
#' @param formula.sigma an object of class formula: a symbolic description of the model to be fitted for the beta dispersion parameter.
#' @param formula.mu an object of class formula: a symbolic description of the model to be fitted for the beta mean parameter.
#' @param formula.nu an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param formula.tau an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param covariates data frame of covariate, Please include intercept as the first column
#' @param rho the list of the coefficient for nu rho[[1]] is the coefficient for x1 and  rho[[2]] is the coefficient for x2.
#' @param delta the list of the coefficient for mu delta[[1]] is the coefficient for x1 and  delta[[2]] is the coefficient for x2.
#' @param kappa the list of the coefficient for sigma kappa[[1]] is the coefficient for x1 and  kappa[[2]] is the coefficient for x2.
#' @param gamma the list of the coefficient for tau gamma[[1]] is the coefficient for x1 and  gamma[[2]] is the coefficient for x2.

#' @export
simulateBivariateZOIBRcov = function(n, theta, formula.mu=~1, formula.sigma=~1, formula.nu=~1, formula.tau=~1, covariates=NULL, rho=NULL, delta=NULL, kappa=NULL, gamma=NULL) {
  # 0. convert covariates to data.matrix so can do %*%
    if(!is.null(covariates)) covariates = data.matrix(covariates)
    
    
    
    # 1. Transformation from rho/delta/kappa to p_i/alpha_i/beta_i
    
    # assign regression parameters their own vectors
    rho1 = rho[[1]]; rho2 = rho[[2]]
    delta1 = delta[[1]]; delta2 = delta[[2]]
    kappa1 = kappa[[1]]; kappa2 = kappa[[2]]
    gamma1 = gamma[[1]]; gamma2 = gamma[[2]]
    
    formula.mu.1 = stats::as.formula(formula.mu[[1]])
    formula.sigma.1 = stats::as.formula(formula.sigma[[1]])
    formula.nu.1 = stats::as.formula(formula.nu[[1]])
    formula.tau.1 = stats::as.formula(formula.tau[[1]])
    
    formula.mu.2 = stats::as.formula(formula.mu[[2]])
    formula.sigma.2 = stats::as.formula(formula.sigma[[2]])
    formula.nu.2 = stats::as.formula(formula.nu[[2]])
    formula.tau.2 = stats::as.formula(formula.tau[[2]])
    
    
    
    terms_mu_1 <- terms(formula.mu.1)
    terms_sigma_1 <- terms(formula.sigma.1)
    terms_nu_1 <- terms(formula.nu.1)
    terms_tau_1 <- terms(formula.tau.1)
    
    terms_mu_2 <- terms(formula.mu.2)
    terms_sigma_2 <- terms(formula.sigma.2)
    terms_nu_2 <- terms(formula.nu.2)
    terms_tau_2 <- terms(formula.tau.2)
    
    #n=length(cell_type)
    
    # Function to construct matrix based on intercept presence
    construct_matrix <- function(terms, covariates, n) {
      if(attr(terms, "intercept") == 1){ # if RHS of y ~ x is == ~ 1 then intercept only model
        if(length(attr(terms, "term.labels"))==0){
          cbind(rep(1, n))
        }else{
          cbind("int"=rep(1, n), covariates[,attr(terms, "term.labels")])
        }
      } else{ # else extract covariates in the model
        covariates[,attr(terms, "term.labels")]
      }
    }
    
    Q1 <- construct_matrix(terms_nu_1, covariates, n)
    Q2 <- construct_matrix(terms_nu_2, covariates, n)
    
    # W (mu)
    W1 <- construct_matrix(terms_mu_1, covariates, n)
    W2 <- construct_matrix(terms_mu_2, covariates, n)
    
    # Z (sigma)
    Z1 <- construct_matrix(terms_sigma_1, covariates, n)
    Z2 <- construct_matrix(terms_sigma_2, covariates, n)
    
    #P (tau)
    P1 <- construct_matrix(terms_tau_1, covariates, n)
    P2 <- construct_matrix(terms_tau_2, covariates, n)
    
    # # Q (nu)
    # Q1 <- construct_matrix(terms_nu, covariates, n)
    # Q2 <- construct_matrix(terms_nu, covariates, n)
    # 
    # # W (mu)
    # W1 <- construct_matrix(terms_mu, covariates, n)
    # W2 <- construct_matrix(terms_mu, covariates, n)
    # 
    # # Z (sigma)
    # Z1 <- construct_matrix(terms_sigma, covariates, n)
    # Z2 <- construct_matrix(terms_sigma, covariates, n)
    # # assign separate design matrices from covariate ls() & add column of 1s for intercept
    # # Q (nu)
    # if(attr(terms_nu, "intercept") == 1){ # if RHS of y ~ x is == ~ 1 then intercept only model
    #   if(length(attr(terms_nu, "term.labels"))==0){
    #     Q1 = cbind(rep(1, n))
    #   }else{
    #     Q1 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_nu, "term.labels"))])
    #   }
    # } else{ # else extract covariates in the model
    #   Q1 = covariates[,all.vars(attr(terms_nu, "term.labels"))]
    # }
    # 
    # if(attr(terms_nu, "intercept") == 1){ # if RHS of y ~ x is == ~ 1 then intercept only model
    #   Q2 = cbind(rep(1, n))
    # } else { # else extract covariates in the model
    #   Q2 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_nu, "term.labels"))])
    # }
    # 
    # # W (mu)
    # if(attr(terms_mu, "intercept") == 1){
    #   W1 = cbind(rep(1, n))
    # } else {
    #   W1 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_mu, "term.labels"))])
    # }
    # 
    # if(attr(terms_mu, "intercept") == 1){
    #   W2 = cbind(rep(1, n))
    # } else {
    #   W2 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_mu, "term.labels"))])
    # }
    # 
    # # Z (sigma)
    # if(attr(terms_sigma, "intercept") == 1){
    #   Z1 = cbind(rep(1, n))
    # } else {
    #   Z1 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_sigma, "term.labels"))])
    # }
    # 
    # if(attr(terms_sigma, "intercept") == 1){
    #   Z2 = cbind(rep(1, n))
    # } else {
    #   Z2 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_sigma, "term.labels"))])
    # }
    # 
    # if(attr(terms_tau, "intercept") == 1){ # if RHS of y ~ x is == ~ 1 then intercept only model
    #   if(length(attr(terms_tau, "term.labels"))==0){
    #     P1 = cbind(rep(1, n))
    #   }else{
    #     P1 = cbind("int"=rep(1, n), covariates[,all.vars(attr(terms_tau, "term.labels"))])
    #   }
    # } else{ # else extract covariates in the model
    #   P1 = covariates[,all.vars(attr(terms_nu, "term.labels"))]
    # }
    
    # ZI-Beta Regression
    # GLM for p
    logit.nu1 = Q1 %*% rho1; logit.nu2 = Q2 %*% rho2
    # GLM for mu
    logit.mu1 = W1 %*% delta1; logit.mu2 = W2 %*% delta2
    # GLM for sigma
    log.sigma1 = Z1 %*% kappa1; log.sigma2 = Z2 %*% kappa2
    # GLM for tau
    log.tau1 = P1 %*% gamma1; log.tau2 = P2 %*% gamma2
    
    # inverse logit function
    inv.logit = function(x) exp(x) / (1 + exp(x))
    
    # nu
    nu1 = exp(logit.nu1); nu2 = exp(logit.nu2)
    # mu
    mu1 = inv.logit(logit.mu1); mu2 = inv.logit(logit.mu2)
    # sigma
    sigma1 = inv.logit(log.sigma1); sigma2 = inv.logit(log.sigma2)
    
    # tau
    tau1 = exp(log.tau1); tau2 = exp(log.tau2)
    
    # shape parameters
    # alpha1 <- mu1*phi1; beta1 <- phi1 - mu1*phi1
    # alpha2 <- mu2*phi2; beta2 <- phi2 - mu2*phi2
    
    
    alpha1 = mu1*(1-sigma1^2)/sigma1^2
    beta1 = (1-mu1)*(1-sigma1^2)/sigma1^2
    p1 = nu1/(1+nu1+tau1)
    q1 = tau1/(1+nu1+tau1)
    
    alpha2 = mu2*(1-sigma2^2)/sigma2^2
    beta2 = (1-mu2)*(1-sigma2^2)/sigma2^2
    p2 = nu2/(1+nu2+tau2)
    q2 = tau2/(1+nu2+tau2)
    
    print(summary(cbind(alpha1,alpha2,beta1,beta2)))
    
    
    checkX0 = FALSE; checkNZPairs = FALSE
    ## check if theta == 0 (independence)
    thetaCheck = theta != 0
    
    # simulate x1 and x2 using ZIB regression margins
    while(!(checkX0 & checkNZPairs)) {
      # simulate u and w, independent uniform random variables
      u = runif(n = n, min = 0, max = 1)                                      # u = F1(x1)
      
      # if theta != 0 simulate v by:
      if(thetaCheck){
        w = runif(n = n, min = 0, max = 1)                                    # w = c1(u,v)
        
        #a = Rmpfr::mpfr(exp(-theta), prec)  
        a = exp(-theta)  # parameterize a = exp(-theta)
        v = as.numeric((-1/theta) * log(1 + (w*(a - 1))/(w + (a^u)*(1 - w)))) # v = F2(x2)
      }
      
      # if theta == 0 simulate v by:
      else {
        v = runif(n = n, min = 0, max = 1)                                    # v = F2(x2)
      }
      
      # **the inverse CDF transformation is done in two stages because doing it in one qbeta gives a warning**
      # if u <= p1, x1u = 0, else transform u into x1u using (u - p1)/(1 - p1)
      # likewise for v and x2v
      x1u = ifelse(u <= p1, 0, (u-p1)/(1-p1-q1))
      x2v = ifelse(v <= p2, 0, (v-p2)/(1-p2-q2))
      x1u = ifelse(u > 1-q1, 1, x1u)
      x2v = ifelse(v > 1-q2, 1, x2v)
      
      # if x1u == 0, x1 = 0, else x1 is the inverse CDF ox x1u
      # likewise for x2v and x2
      x1 = ifelse((x1u == 0)|(x1u == 1), x1u, qbeta(x1u, alpha1, beta1))      # use the inverse CDF to solve for x1
      x2 = ifelse((x2v == 0)|(x2v == 1), x2v, qbeta(x2v, alpha2, beta2))      # use the inverse CDF to solve for x2
      
      # check0 = length(which(x1 != 0)) < 3 | length(which(x2 != 0)) < 3 # each x needs >= 3 non-zero obs to estimate parameters
      checkX0 = length(which(x1 != 0)) > 3 | length(which(x2 != 0)) > 3 # each x needs >= 3 non-zero obs to estimate parameters
      checkNZPairs = length(which(x1!=0 & x2!=0)) > 1                   # at least 2 non-zero pairs are needed to stable est
    }
    
    return(x = data.frame(x1, x2))

 
}

