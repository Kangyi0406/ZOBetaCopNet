#' The Zero-One-Inflated Beta Distribution
#'
#' Density function for the zero-one-inflated beta distribution with zero-one-inflation probability equal to \code{p}, one-inflation probability equal to \code{q}, beta distribution alpha equal to \code{alpha} and beta distribution beta equal to \code{beta}.
#'
#' @param x vector of quantiles.
#' @param nu vector of zero-inflation probabilities.
#' @param tau vector of one-inflation probabilities.
#' @param alpha vector of alpha in beta distribution
#' @param beta vector of beta in beta distribution
#'
#' @return \code{dzib} returns the density.
#'
#' @export
#'
#' @examples
#' dzib(x = 0.8, p = 0.2, q=0.3, alpha = 2, beta = 3)
dzib = function(x, p, q, alpha, beta){
  if(class(x)!="numeric" | all(x >= 0) == FALSE | all(x <= 1) == FALSE) {
    stop("ERROR: x must be a numeric vector of zero-inflated beta random variables with range [0,1].")
  } else if(class(p) !="numeric" | all(p >= 0) == FALSE | all(p < 1) == FALSE) {
    stop("ERROR: p must be a numeric vector of zero-inflation probabilities with range [0,1).")
  } else if(class(q) !="numeric" | all(q >= 0) == FALSE | all(q < 1) == FALSE) {
    stop("ERROR: q must be a numeric vector of one-inflation probabilities with range [0,1).")
  } else if(class(alpha1) !="numeric"  |
            class(alpha2) !="numeric" ) {
    stop("ERROR: alpha1 and alpha2 must be numeric.")
  } else if(class(beta1) !="numeric"  |
            class(beta2) !="numeric" ) {    stop("ERROR: beta1 and beta2 must be numeric.")  } 
  get_fx = function(i){
    if(x[i]==0){
      return(p[i])
    }else if(x[i]==1){
      return(q[i])
    }else{
      return((1 - p[i] - q[i]) * stats::dbeta(x[i], alpha[i], beta[i]))
    }
  }
  fx = sapply(1:length(x),get_fx)
  
  
  return(fx)
}
#' 