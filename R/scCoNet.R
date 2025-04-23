#' Estimate Marginal Parameters and Network Connections for Large-Scale Data
#'
#' This function computes the maximum likelihood estimates (MLEs) of marginal parameters for each feature in a large dataset, optionally adjusting for covariates. It supports parallel computation and stores the estimated parameters. Subsequently, pairwise associations between features are quantified using likelihood ratio tests, returning both the test statistics and significance indicators.
#'
#' @param data A data frame where rows represent samples and columns represent features (e.g., genes).
#' @param covars A data frame of covariates used to adjust the marginal models. Row order must match that of \code{data}, and row names should correspond to sample identifiers.
#' @param ncores Integer. Number of cores to use for parallel computation. Default is 1 (no parallelization).
#' @param formula.mu A formula specifying the model for the beta distribution's mean (\eqn{\mu}) parameter. Default is \code{y ~ 1}, indicating no covariates.
#' @param formula.sigma A formula specifying the model for the beta distribution's dispersion (\eqn{\sigma}) parameter. Default is \code{y ~ 1}.
#' @param formula.nu A formula specifying the model for the zero-one inflation parameter (\eqn{\nu}). Default is \code{y ~ 1}.
#' @param formula.tau A formula specifying the model for the complementary zero-one inflation parameter (\eqn{\tau}). Default is \code{y ~ 1}.
#' @param sig_level Numeric. Significance threshold (after Benjamini–Hochberg adjustment) used to determine which associations are retained. Default is 0.05.
#' @param sig_connect Logical. If \code{TRUE}, only statistically significant associations are returned; otherwise, all estimated connections are included. Default is \code{FALSE}.
#' @import progress
#' @import doParallel
#' @import foreach
#' @import progress
#' @return A list containing two elements:
#'   \code{marginal_parameter}{A list of estimated marginal parameters for each feature. Each entry is a matrix with columns \code{q}, \code{p}, \code{alpha}, and \code{beta}.
#'   \code{data_frame}{A data frame summarizing two-stage estimated theta and the testing. Columns include: feature indices, estimated theta, transformed Spearman correlation, likelihood ratio statistic, raw p-value, BH-adjusted p-value, and significance indicator.
#' @examples
#' data(gene_data)
#' data = gene_data$gen_data
#' cov = gene_data$cov
#' cov$name_data = factor(cov$name_data)
#' cov$year = as.numeric(cov$year)
#' result = scCoNet(data,covars= data.frame(cov[,2:3]), ncores = 1, formula.mu = y~name_data+year,formula.sigma = y~name_data+year, formula.nu = y~name_data+year)
#' @export

scCoNet = function(data, covars=NULL, ncores = 1, sig_level = 0.05, sig_connect = FALSE, formula.mu = y~1, formula.sigma = ~1, formula.nu = ~ 1, formula.tau = ~ 1){
  data =  as.data.frame(data)
  results = list()
  if (is.null(covars)){
    get_marginal <- function(x){
      df <- data.frame(y = stats::na.omit(as.numeric(x)))
      temp <- marginalMLE(formula.mu = y ~ 1, formula.sigma = ~ 1, formula.nu = ~ 1, formula.tau = ~ 1, df = df)
      return(temp$fitted)

    }
    for(i in 1:dim(data)[2]){
      results[[i]] = get_marginal(data[,i])
      #print(i)
    }
    cat("All columns of Marginal MLE processed.\n")

  }else{
    get_marginal= function(x){
      x =as.numeric(x)
      df = cbind(x, covars)
      colnames(df) = c("y",colnames(covars))
      df = stats::na.omit(df)
      temp = marginalMLE(formula.mu = formula.mu, formula.sigma =formula.sigma, formula.nu = formula.nu, formula.tau = formula.tau, df = df)
      #temp = marginalMLE(formula.mu = y ~ source, formula.phi = ~ source, formula.p = ~ source, formula.q = ~ source, df = df)
      return(temp$fitted)
    }
    for(i in 1:dim(data)[2]){
      results[[i]] = get_marginal(data[,i])
      #print(i)
    }
    cat("All columns of Marginal MLE processed.\n")


  }
  frank_theta_to_spearman <- function(theta) {
    frankCopula <- copula::archmCopula(family = "frank", param = theta)
    spearman_rho <- copula::rho(frankCopula)
    return(spearman_rho)
  }
  compute_theta <- function(pair,results,data){
    col1 <- as.numeric(data[,pair[1]])
    col2 <- as.numeric(data[,pair[2]])
    # col1 <- as.numeric(data[,2])
    # col2 <- as.numeric(data[,4])


    x <- stats::na.omit(as.data.frame(cbind(col1, col2)))



    p1 = results[[pair[1]]]$p
    q1 = results[[pair[1]]]$q
    alpha1 = results[[pair[1]]]$alpha
    beta1 = results[[pair[1]]]$beta
    p2 = results[[pair[2]]]$p
    q2 = results[[pair[2]]]$q
    alpha2 = results[[pair[2]]]$alpha
    beta2 = results[[pair[2]]]$beta
    theta_get = thetaTSMLE(x,
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
    theta_test = lrtTheta(x,
                          p1 = results[[pair[1]]]$p,
                          q1 = results[[pair[1]]]$q,
                          alpha1 = results[[pair[1]]]$alpha,
                          beta1 = results[[pair[1]]]$beta,
                          p2 = results[[pair[2]]]$p,
                          q2 = results[[pair[2]]]$q,
                          alpha2 = results[[pair[2]]]$alpha,
                          beta2 = results[[pair[2]]]$beta,
                          theta = theta_get
    )
    # spearman_test = cor.test(x[,1],x[,2],method = "spearman",use="pairwise.complete.obs")
    theta_est = frank_theta_to_spearman(theta_get)
    if(theta_test[1] <=0){
      theta_est = 0
    }
    return(c(pair[1], pair[2], theta_get, theta_est, theta_test))
  }
  num_cols = ncol(data)
  # Create all combinations of column pairs
  all_pairs <- utils::combn(1:num_cols, 2, simplify = FALSE)
  # Process each chunk

  # Set up cluster for this chunk
  if(ncores==1){
    results_theta = data.frame(do.call(rbind,lapply( all_pairs, compute_theta, results =  results, data = data)))
    colnames(results_theta) = c("Column1", "Column2", "Copula","Copula_Spearman", "Test_Statistics","Copula_P_value")
    results_theta = data.frame(results_theta)
    results_theta$Adjusted_P_copula <- stats::p.adjust(results_theta$Copula_P_value, method = "BH")
    results_theta$Significant_copula <- results_theta$Adjusted_P_copula < sig_level

    if(sig_connect){
      results_theta = results_theta[results_theta$Significant_copula ==1,]
    }
    return(list(marginal_parameter = results, data_frame=results_theta))

  }


  ###multi-cores situation

  cl <- parallel::makeCluster(ncores, type = "SOCK")
  doSNOW::registerDoSNOW(cl)
  # Export required data and objects to all workers
  parallel::clusterExport(cl, varlist = c("data", "results"), envir = environment())
  # Load your package on all workers:
  parallel::clusterCall(cl, function() library(scCoNet))

  # Progress bar for this chunk
  chunk_pb <- progress::progress_bar$new(
    format = "Processing :current/:total [:bar] :percent eta: :eta",
    total = length(all_pairs),
    clear = FALSE,
    width = 60
  )

  progress <- function(n) chunk_pb$tick()
  opts <- list(progress = progress)


  results_theta <- foreach::foreach(pair = all_pairs,
                           .combine = rbind,
                           .packages = c("copula","data.table"),
                           .export = c("thetaTSMLE","lrtTheta","bivariateDensityZOIB","pzib","dzib","frankConditionalV","frankConditionalU","frankLogLik","frankDensity","frankCopula"),
                           .options.snow = opts,
                           .errorhandling = 'pass',
                           .inorder = T,
                           .multicombine = TRUE) %dopar% {
                             compute_theta(pair, results =  results, data = data)
                           }

  colnames(results_theta) = c("Column1", "Column2", "Copula","Copula_Spearman", "Test_Statistics","Copula_P_value")
  results_theta = data.frame(results_theta)
  results_theta$Adjusted_P_copula <- stats::p.adjust(results_theta$Copula_P_value, method = "BH")
  results_theta$Significant_copula <- results_theta$Adjusted_P_copula < sig_level

  if(sig_connect){
    results_theta = results_theta[results_theta$Significant_copula ==1,]
  }

    # Stop the cluster
  return(list(marginal_parameter = results, data_frame=results_theta))
  parallel::stopCluster(cl)



}
