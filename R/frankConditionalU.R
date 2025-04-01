#' The Frank Conditional Copula
#'
#' Frank conditional copula function.
#'
#' @param u vector of uniform random variables.
#' @param v vector of uniform random variables.
#' @param theta copula dependence parameter.

#'
#' @return \code{frankConditionalU} returns the frank conditional copula distribution function, conditional on u.
#'
#' @export
#'
#' @examples
#' frankConditionalU(u = 0.35, v = 0.5, theta = 2)
#' 
frankConditionalU = function(u, v, theta){
  if(class(u)!="numeric" | all(u >= 0) == FALSE | all(u <= 1) == FALSE) {
    stop("ERROR: u must be a numeric vector of Uniform(0,1) random variables.")
  } else if(class(v) !="numeric" | all(v >= 0) == FALSE | all(v <= 1) == FALSE) {
    stop("ERROR: v must be a numeric vector of Uniform(0,1) random variables.")
  }
  a = exp(-theta)
  
  
  C1uv = as.numeric((a^u * (a^v - 1)) / ((a^u - 1) * (a^v - 1) + (a - 1)))
  
  return(C1uv)
}