#' The Zero-One-Inflated Beta Distribution
#'
#' Density function for the zero-one-inflated beta distribution with zero-one-inflation probability equal to \code{p}, one-inflation probability equal to \code{q}, beta distribution alpha equal to \code{alpha} and beta distribution beta equal to \code{beta}.
#'
#' @param x vector of quantiles.
#' @param p vector of zero-inflation probabilities.
#' @param q vector of one-inflation probabilities.
#' @param alpha vector of alpha in beta distribution
#' @param beta vector of beta in beta distribution
#'
#' @return \code{dzib} returns the density.
#'
#'
#' @export

dzib = function(x, p, q, alpha, beta){
  if (!is.numeric(x) || !all(x >= 0) || !all(x <= 1)) {
    stop("ERROR: x must be a numeric vector of zero-inflated beta random variables with range [0,1].")
  } else if (!is.numeric(p) || !all(p >= 0) || !all(p < 1)) {
    stop("ERROR: p must be a numeric vector of zero-inflation probabilities with range [0,1).")
  } else if (!is.numeric(q) || !all(q >= 0) || !all(q < 1)) {
    stop("ERROR: q must be a numeric vector of one-inflation probabilities with range [0,1).")
  } else if (!is.numeric(alpha)) {
    stop("ERROR: alpha must be numeric.")
  } else if (!is.numeric(beta)) {
    stop("ERROR: beta must be numeric.")
  }

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
