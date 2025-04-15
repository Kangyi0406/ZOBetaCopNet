#' Bivariate Zero-One-Inflated Beta Density
#'
#' Bivariate frank copula density function with zero-inflated beta margins.
#'
#' @param x1 first vector of relative abundances in the (x1, x2) pair.
#' @param p1 zero-inflation probability for x1.
#' @param q1 one-inflation probability for x1.
#' @param alpha1 alpha parameter in beta distribution for x1.
#' @param beta1 beta parameter in beta distribution for x1.
#' @param x2 second vector of relative abundances in the (x1, x2) pair.
#' @param p2 zero-inflation probability for x2.
#' @param q2 one-inflation probability for x2.
#' @param alpha2 alpha parameter in beta distribution for x2.
#' @param beta2 beta parameter in beta distribution for x2.
#' @param theta copula dependence parameter.
#'
#' @return \code{bivariateDensityZOIB} returns the bivariate zero-one-inflated beta density function using a frank copula.
#'
#' @export

bivariateDensityZOIB = function(x1, p1, q1, alpha1, beta1, x2, p2, q2, alpha2, beta2, theta){
  if (!is.numeric(x1) || !all(x1 >= 0) || !all(x1 <= 1) ||
      !is.numeric(x2) || !all(x2 >= 0) || !all(x2 <= 1)) {
    stop("ERROR: x1 and x2 must be numeric vectors of zero-inflated beta random variables with range [0,1].")
  } else if (!is.numeric(p1) || !all(p1 >= 0) || !all(p1 < 1) ||
             !is.numeric(p2) || !all(p2 >= 0) || !all(p2 < 1)) {
    stop("ERROR: p1 and p2 must be numeric vectors of zero-inflation probabilities with range [0,1).")
  } else if (!is.numeric(q1) || !all(q1 >= 0) || !all(q1 < 1) ||
             !is.numeric(q2) || !all(q2 >= 0) || !all(q2 < 1)) {
    stop("ERROR: q1 and q2 must be numeric vectors of one-inflation probabilities with range [0,1).")
  } else if (!is.numeric(alpha1) || !is.numeric(alpha2)) {
    stop("ERROR: alpha1 and alpha2 must be numeric.")
  } else if (!is.numeric(beta1) || !is.numeric(beta2)) {
    stop("ERROR: beta1 and beta2 must be numeric.")
  } else if (!is.numeric(theta)) {
    stop("ERROR: theta must be a numeric variable.")
  }



  fx1 = dzib(x1, p1, q1, alpha1, beta1) # density of x1
  fx2 = dzib(x2, p2, q2, alpha2, beta2) # density of x2

  # cdf
  Fx1 = pzib(x1, p1, q1, alpha1, beta1) # cdf at x1
  Fx2 = pzib(x2, p2, q2, alpha2, beta2) # cdf at x2

  # empty joint distribution vector
  fx1x2 = vector(length = length(x1))

  #the function return in continous set or not
  notinZeroOne = function(x){
    return(x != 0 & x != 1)
  }

  # subset obs. by each combo
  s11 = which(notinZeroOne(x1) & notinZeroOne(x2)) # number of x1 not in {0,1} & x2 not in {0,1}
  s01 = which(x1 == 0 & notinZeroOne(x2)) # number of x1 == 0 & x2 not in {0,1}
  s10 = which(notinZeroOne(x1) & x2 == 0) # number of x1 not in {0,1} & x2 == 0
  s21 = which(x1 == 1 & notinZeroOne(x2)) # number of x1 == 1 & x2 not in {0,1}
  s12 = which(notinZeroOne(x1) & x2 == 1) # number of x1 not in {0,1} & x2 == 1
  s00 = which(x1 == 0 & x2 == 0) # number of x1 == 0 & x2 == 0
  s22 = which(x1 == 1 & x2 == 1) # number of x1 == 1 & x2 == 1
  s20 = which(x1 == 1 & x2 == 0) # number of x1 == 1 & x2 == 0
  s02 = which(x1 == 0 & x2 == 1) # number of x1 == 0 & x2 == 1



  # joint distribution -- using frank copula
  fx1x2[s11] = frankDensity(Fx1[s11], Fx2[s11], theta)*fx1[s11]*fx2[s11] # joint density for s1
  fx1x2[s01] = fx2[s01]*frankConditionalV(p1[s01], Fx2[s01], theta)    # joint density for s2
  fx1x2[s10] = fx1[s10]*frankConditionalU(Fx1[s10], p2[s10], theta)    # joint density for s3
  fx1x2[s21] = fx2[s21]*(frankConditionalV(rep(1,length(s21)), Fx2[s21], theta) - frankConditionalV(1-q1[s21], Fx2[s21], theta))    # joint density for s4
  fx1x2[s12] = fx1[s12]*(frankConditionalU(Fx1[s12], rep(1,length(s12)), theta) - frankConditionalU(Fx1[s12], 1-q2[s12], theta))    # joint density for s5
  fx1x2[s00] = frankCopula(Fx1[s00], Fx2[s00], theta)                  # joint density for s6
  fx1x2[s22] = frankCopula(rep(1,length(s22)), rep(1,length(s22)), theta) - frankCopula(rep(1,length(s22)), 1-q2[s22], theta) - frankCopula(1-q1[s22], 1, theta) + frankCopula(1-q1[s22], 1-q2[s22], theta)                  # joint density for s7
  fx1x2[s20] = frankCopula(rep(1,length(s20)), p2[s20], theta)- frankCopula(1-q1[s20], p2[s20], theta)                  # joint density for s8
  fx1x2[s02] = frankCopula(p1[s02], rep(1,length(s02)), theta) - frankCopula(p1[s02], 1-q2[s02], theta)
  # joint density for s9


  return(fx1x2)
}


