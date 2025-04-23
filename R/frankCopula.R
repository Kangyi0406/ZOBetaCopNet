#' The Frank Copula
#'
#' Frank copula function.
#'
#' @param u vector of uniform random variables.
#' @param v vector of uniform random variables.
#' @param theta copula dependence parameter.
#'
#' @return \code{frankCopula} returns the frank copula distribution function.
#'
#'
#' @export
#'
frankCopula = function(u, v, theta){
  if (!is.numeric(u) || !all(u >= 0) || !all(u <= 1)) {
    stop("ERROR: u must be a numeric vector of Uniform(0,1) random variables.")
  } else if (!is.numeric(v) || !all(v >= 0) || !all(v <= 1)) {
    stop("ERROR: v must be a numeric vector of Uniform(0,1) random variables.")
  } else if (!is.numeric(theta)) {
    stop("ERROR: theta must be a numeric variable.")
  }
  a = exp(-theta)

  Cuv = (-1/theta)*as.numeric(log(1 + (((a^u - 1) * (a^v - 1))/(a - 1))))

  return(Cuv)
}
