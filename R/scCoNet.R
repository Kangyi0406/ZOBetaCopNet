
#' Large Data get the MLE of marginal parameter and get the estimate
#'
#' Get the MLE of marginal parameter for the Large Dataset and save it and then estimate the parameter
#'
#'
#' @param abd data frame of microbial relative abundances. Rows are samples and columns are taxa.
#' @param covars data frame of covariates to adjust for in the marginal regression models. Rows should be in the same order as \code{abd} and rownames should match.
#' @param ncores number of cores for parallelization. Default value is 1 (no parallel computing).
#' @param wd_Marginal working directory file name for saving the rds file related to marginal MLE parameter, the default is NULL
#' @param wd_estimate working directory file name for saving the rds file related to copula estimator and p-value,, the default is NULL
#' @param formula.sigma an object of class formula: a symbolic description of the model to be fitted for the beta dispersion parameter.
#' @param formula.mu an object of class formula: a symbolic description of the model to be fitted for the beta mean parameter.
#' @param formula.nu an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param formula.tau an object of class formula: a symbolic description of the model to be fitted for the zero-one inflation probability parameter.
#' @param sig_level the significant level for select the connection based on the B-H adjusted p-value, defualt is 0.05.
#' @param sig_connect the argument for whether only return the significant connection or not. The default value is FALSE, it return all the connections.
#'
#' @return \code{scCoNet} returns dataframe (Column 1 index, Column 2 index, estimated theta, transformed spearman, LR statistic, p-value, p-value after B-H adjustment)
#' @import progress
#' @import doParallel
#' @import foreach
#' @import progress
#' @examples
#' cell_type = "Inh
#' data <- as.data.frame(readRDS(paste0("/storage/Network_data/data_combine/beta_protein/beta_",cell_type,"_use.rds")))[,1:10]
#' cov = as.data.frame(readRDS(paste0("/storage/Network_data/data_combine/cell_cov/cov_",cell_type,".rds")))
#' cov$name_data = factor(cov$name_data)
#' cov$year = as.numeric(cov$year)
#' scCoNet(data[,1:5],covars= data.frame(cov[,2:3]), ncores = 2, formula.mu = y~name_data+year)
#' @export

scCoNet = function(abd, covars=NULL, ncores = 1, sig_level = 0.05, sig_connect = FALSE, wd_Marginal=NULL, wd_estimate= NULL, formula.mu = y~1, formula.sigma = ~1, formula.nu = ~ 1, formula.tau = ~ 1){
  data =  as.data.frame(abd)
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
    if(!is.null(wd_Marginal)){
      saveRDS(results, wd_Marginal)
      cat("All columns of Marginal MLE processed and results saved.\n")
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
    if(!is.null(wd_Marginal)){
      saveRDS(results, wd_Marginal)
      cat("All columns of Marginal MLE processed and results saved.\n")
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
  print(results_theta)

  colnames(results_theta) = c("Column1", "Column2", "Copula","Copula_Spearman", "Test_Statistics","Copula_P_value")
  results_theta = data.frame(results_theta)
  results_theta$Adjusted_P_copula <- stats::p.adjust(results_theta$Copula_P_value, method = "BH")
  results_theta$Significant_copula <- results_theta$Adjusted_P_copula < sig_level

  if(sig_connect){
    results_theta = results_theta[results_theta$Significant_copula ==1,]
  }

    # Stop the cluster
    # Save results for this chunk
    if(is.null(wd_estimate)){
      return(results_theta)
    }
    saveRDS(results_theta, wd_estimate)
    parallel::stopCluster(cl)

    message("All estimations processed and saved.\n")


}
