#' copnet
#'
#'  Implementation of Copulas with Mixture Margins for Covariation Networks using two-stage estimation and testing.
#'
#' @param abd data frame of microbial relative abundances. Rows are samples and columns are taxa.
#' @param covars data frame of covariates to adjust for in the marginal regression models. Rows should be in the same order as \code{abd} and rownames should match.
#' @param ncores number of cores for parallelization. Default value is 1 (no parallel computing).
#'
#' @return \code{copnet} returns the two-stage estimate of \eqn{\theta} the LR statistics and the p-value.
#' @export
#'
copnet = function(abd, covars = NULL, ncores = 1) {
  #covars = data.frame(covariates)


  if(!("data.frame"%in% class(abd))) {
    stop("ERROR: abd must be a data frame of relative abundances.")
  } else if(!is.null(covars) & !("data.frame" %in% class(covars))) {
    stop("ERROR: covars must be a data frame of covariates covariates for the marginal regression models.")
  }


  if(any(abd < 0) | any(abd > 1)) {
    stop("ERROR: the entries of abd must be relative abundances between [0,1].")
  } else if(is.null(covars) == FALSE ) {
    if(all.equal(rownames(abd), rownames(covars)) != TRUE) {
      stop("ERROR: the rownames of abd and covars must match and be in the same order.")
    }
  }


  if(ncores > 1){
    message("Setting up paralleization")
    future::plan(future::multisession, workers = ncores)
  }

  #message("Fitting marginal models.")
  if(is.null(covars)){
    mMLE = furrr::future_map(colnames(abd),
                             function(y, abd_mat) {
                               df = data.frame(y = abd_mat[,y])
                               out = marginalMLE(formula.mu = y ~ 1, formula.sigma = ~ 1, formula.nu = ~ 1,formula.tau = ~ ., df = df)
                               return(out)
                             },
                             abd_mat = abd,
                             .progress = TRUE, .options = furrr::furrr_options(seed = T) )
  } else {
    mMLE = furrr::future_map(colnames(abd),
                             function(y, abd_mat, vars) {
                               df = cbind(y = abd_mat[,y], vars)
                               out = marginalMLE(formula.mu = y ~ ., formula.sigma = ~ ., formula.nu = ~  ., formula.tau = ~ ., df = df)
                               return(out)
                             },
                             abd_mat = abd, vars = covars,
                             .progress = TRUE, .options = furrr::furrr_options(seed = T) )
  }


  message("Estimating copula dependence parameters.")
  p1 = mMLE[[1]]$fitted$p
  p2 = mMLE[[2]]$fitted$p
  q1 = mMLE[[1]]$fitted$q
  q2 = mMLE[[2]]$fitted$q
  alpha1 = mMLE[[1]]$fitted$alpha
  alpha2 = mMLE[[2]]$fitted$alpha
  beta1 = mMLE[[1]]$fitted$beta
  beta2 = mMLE[[2]]$fitted$beta


  theta = thetaTSMLE(x = abd,
                     lower = -30,
                     upper = 30,
                     p1 = p1,
                     q1 = q1,
                     alpha1 = alpha1,
                     beta1 = beta1,
                     p2 = p2,
                     q2 = q2,
                     alpha2 = alpha2,
                     beta2 = beta2

  )


    message("Performing likelihood ratio tests.")
    lR.ts = lrtTheta(response = abd,
                     theta = theta,
                     p1 = p1,
                     q1 = q1,
                     alpha1 = alpha1,
                     beta1 = beta1,
                     p2 = p2,
                     q2 = q2,
                     alpha2 = alpha2,
                     beta2 = beta2,
                     thetaNull = 0)


  return(c(theta,lR.ts))
}

