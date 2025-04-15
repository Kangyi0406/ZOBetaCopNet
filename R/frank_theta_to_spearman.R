#' Convert theta in Frank copula function to spearman correlation
#'
#' Convert theta in Frank copula function to spearman correlation
#' @param theta The copula parameter in Frank copula function
#'
#' @return \code{frank_theta_to_spearman} returns the corresponding spearman correlation for copula parameter in Frank copula function
#'
#' @export

frank_theta_to_spearman <- function(theta) {
  frankCopula <- copula::archmCopula(family = "frank", param = theta)
  spearman_rho <- copula::rho(frankCopula)
  return(spearman_rho)
}
