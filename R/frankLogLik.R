#' Log-likelihood for the Bivariate Zero-One-Inflated Beta Density
#'
#' Log-likelihood for the bivariate frank copula density function with a zero-one-inflated beta margins.
#'
#'
#' @param x data frame with two columns for bivariate pair (x1, x2).
#' @param p1 zero-inflation probability for x1.
#' @param q1 one-inflation probability for x1.
#' @param alpha1 alpha parameter in beta distribution for x1.
#' @param beta1 beta parameter in beta distribution for x1.
#' @param p2 zero-inflation probability for x2.
#' @param q2 one-inflation probability for x2.
#' @param alpha2 alpha parameter in beta distribution for x2.
#' @param beta2 beta parameter in beta distribution for x2.
#' @param theta copula dependence parameter.
#'
#' @export

frankLogLik <- function(theta, p1, q1, alpha1, beta1, p2, q2, alpha2, beta2, x) {
  if(class(theta) !="numeric") {
    stop("ERROR: theta must be a numeric variable.")
  } else if(class(p1) !="numeric" | all(p1 >= 0) == FALSE | all(p1 < 1) == FALSE |
            class(p2) !="numeric" | all(p2 >= 0) == FALSE | all(p2 < 1) == FALSE) {
    stop("ERROR: p1 and p2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if(class(q1) !="numeric" | all(q1 >= 0) == FALSE | all(q1 < 1) == FALSE |
            class(q2) !="numeric" | all(q2 >= 0) == FALSE | all(q2 < 1) == FALSE) {
    stop("ERROR: q1 and q2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if(class(alpha1) !="numeric"  |
                class(alpha2) !="numeric" ) {
    stop("ERROR: alpha1 and alpha2 must be numeric.")
  } else if(class(beta1) !="numeric"  |
            class(beta2) !="numeric" ) {
    stop("ERROR: beta1 and beta2 must be numeric.")  
    } else if(!("data.frame" %in% class(x)) | ncol(x) != 2) {
    stop("ERROR: x must be a data frame with two columns containing the bivariate relative abundances.")
  }
  
  # parameterize a = exp(-theta)
  a = exp(-theta)
  # define x1 and x2
  x1 <- x[,1]; x2 <- x[,2]
  
  # Check if theta != 0
  thetaCheck <- theta != 0
  
  if(thetaCheck){ # theta != 0, use copula log-likelihood
    ll = sum(log(bivariateDensityZOIB(x1, p1, q1, alpha1, beta1, x2, p2, q2, alpha2, beta2, theta, prec = 200)))
  }else { 
    # theta == 0, model reduce to independence-- product of marginal pdfs
    # PDF for each x_i
    f1 = dzib(x1, p1, q1, alpha1, beta1); f2 = dzib(x2, p2, q2, alpha2, beta2)
    
    ll = sum(log(f1*f2))
  }
  
  return(as.numeric(ll))
}