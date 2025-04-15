#' Two-stage Likelihood Ratio Test
#'
#' Two-stage likelihood ratio test for \eqn{\theta}.
#'
#' @param response data frame with two columns for bivariate pair (x1, x2).
#' @param p1 zero-inflation probability for x1.
#' @param q1 one-inflation probability for x1.
#' @param alpha1 alpha parameter in beta distribution for x1.
#' @param beta1 beta parameter in beta distribution for x1.
#' @param p2 zero-inflation probability for x2.
#' @param q2 one-inflation probability for x2.
#' @param alpha2 alpha parameter in beta distribution for x2.
#' @param beta2 beta parameter in beta distribution for x2.
#' @param theta copula dependence parameter.
#' @param thetaVar jackknife variance of \eqn{\theta}.
#' @param thetaNull value of \eqn{\theta} under the null hypothesis. Must be 0 for this implementation.
#'
#' @return \code{lrtTheta} returns the two-stage likelihood ratio test statistic and the corresponding p-value.
#'
#' @export
#'
lrtTheta <- function(response, p1, q1, alpha1, beta1, p2, q2, alpha2, beta2, theta, thetaVar = NULL, thetaNull = 0){
  # calculate the log-likelihood under the null (theta = 0)
  loglik.null <- frankLogLik(theta = thetaNull,
                             p1 = p1, q1 = q1, alpha1 = alpha1, beta1 = beta1,
                             p2 = p2, q2 = q2, alpha2 = alpha2, beta2 = beta2,
                             x = response)

  # calculate the log-likelihood under the alternative (theta = MLE)
  loglik.alt <- frankLogLik(theta = theta,
                            p1 = p1, q1 = q1, alpha1 = alpha1, beta1 = beta1,
                            p2 = p2, q2 = q2, alpha2 = alpha2, beta2 = beta2,
                            x = response)

  # IFF thetaNull != 0: calculate scaling factor
  if(thetaNull != 0){
    # sample size
    n = nrow(response)

    jkV.theta = n*thetaVar

    ## approximate two-stage information for theta
    Idd = (-1/n)*maxLik::hessian(maxLik::maxLik(logLik = frankLogLik,
                                                start = c(theta = theta), x = response,
                                                p1 = p1, q1 = q1, alpha1 = alpha1, beta1 = beta1,
                                                p2 = p2, q2 = q2, alpha2 = alpha2, beta2 = beta2))

    # LRT test statistic
    lr.ts <- -2*(loglik.null - loglik.alt)*(jkV.theta*Idd)^-1
    #print(loglik.null)
  }

  # else (thetaNull == 0): scaling factor reduces to 1
  else lr.ts <- -2*(loglik.null - loglik.alt)
  p_value <- stats::pchisq(lr.ts, df = 1, lower.tail = FALSE)


  return(c(as.numeric(lr.ts),p_value))
  #return(as.numeric(lr.ts))
}
