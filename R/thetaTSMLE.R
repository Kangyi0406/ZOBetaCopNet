#' Two-stage Estimation of \eqn{\theta}
#'
#' Two-stage estimation of copula dependence parameter \eqn{\theta}.
#'
#' @param x data frame with two columns for bivariate pair (x1, x2).
#' @param lower the lower end point of the interval to be searched in univariate optimizer.
#' @param upper the lower end point of the interval to be searched in univariate optimizer.
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

#' @return \code{thetaTSMLE} returns the two-stage estimate of \eqn{\theta}.
#'
#' @export
#' 
thetaTSMLE <- function(x, lower, upper, p1, q1, alpha1, beta1, p2, q2, alpha2, beta2) {
  # x = abd
  # lower = -30
  # upper = 30
  # p1 = p1
  # q1 = q1
  # alpha1 = alpha1
  # beta1 = beta1
  # p2 = p2
  # q2 = q2
  # alpha2 = alpha2
  # beta2 = beta2
  x1 <- x[,1]; x2 <- x[,2]

  f <- function(theta){
    density_need = bivariateDensityZOIB(x1, p1, q1, alpha1, beta1, x2, p2, q2, alpha2, beta2, theta)
    #The threshold for ignore the unstable value
    threshold <- 1e-6
    density_final =  density_need[density_need>threshold]
    ll <- -sum(log(density_final))
    #print(ll)
  }
  
  # Numerical analysis to find MLE of theta
  #thetaTSML <- optimize(f, c(lower, upper), maximum = TRUE)[1]
  thetaTSML <- optim(par = 0, fn = f, 
                     method = "Brent", lower = lower, upper = upper)[1]
  
  
  return(as.numeric(thetaTSML))
}
