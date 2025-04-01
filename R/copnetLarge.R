
#' Large Data get the MLE of marginal parameter and get the estimate
#'
#' Get the MLE of marginal parameter for the Large Dataset and save it and then estimate the parameter
#'
#'
#' @param abd data frame of microbial relative abundances. Rows are samples and columns are taxa.
#' @param covars data frame of covariates to adjust for in the marginal regression models. Rows should be in the same order as \code{abd} and rownames should match.
#' @param covar_formula_list list of formula.mu, formula.sigma, formula.nu, formula.tau
#' @param ncores number of cores for parallelization. Default value is 1 (no parallel computing).
#' @param wd_Marginal working directory file name for saving the rds file related to marginal MLE parameter
#' @param wd_estimate working directory file name for saving the rds file related to copula estimator and p-value
#' 
#' 
#'
#' @return \code{copnetLarge} returns (column numbers in the pair, estimated theta, LR statistic, p-value)
#'
#' @export
#' 

copnetLarge = function(abd, covars=NULL, ncores = 1, wd_Marginal, wd_estimate){
  data =  as.data.frame(abd)
  if (covars == NULL){
    result = list()
    get_marginal <- function(x){
      df <- data.frame(y = na.omit(as.numeric(x)))
      temp <- marginalMLE(formula.mu = y ~ 1, formula.sigma = ~ 1, formula.nu = ~ 1, formula.tau = ~ 1, df = df)
      return(temp$fitted)
      
    }
    for(i in 1:dim(data)[2]){
      result[[i]] = get_marginal(data[,i])
      #print(i)
    }
    saveRDS(result, wd_Marginal)
    cat("All columns of Marginal MLE processed and results saved.\n")
  }else{
    get_marginal= function(x){
      x =as.numeric(x)
      df = cbind(x, covars)
      colnames(df) = c("y",colnames(covars))
      temp = marginalMLE(formula.mu = covar_formula_list$mu, formula.sigma =covar_formula_list$sigma, formula.nu = covar_formula_list$mu$mu, formula.tau = covar_formula_list$tau, df = df)
      #temp = marginalMLE(formula.mu = y ~ source, formula.phi = ~ source, formula.p = ~ source, formula.q = ~ source, df = df)
      return(temp$fitted)
    }
    for(i in 1:dim(data)[2]){
      result[[i]] = get_marginal(data[,i])
      #print(i)
    }
    saveRDS(result, wd_Marginal)
    cat("All columns of Marginal MLE processed and results saved.\n")
    
  
  }
  
  # get the marginal MLE
  results = readRDS( wd_Marginal)
  compute_theta <- function(pair){
    col1 <- as.numeric(data[,pair[1]])
    col2 <- as.numeric(data[,pair[2]])
    # col1 <- as.numeric(data[,2])
    # col2 <- as.numeric(data[,4])
    
    
    x <- na.omit(as.data.frame(cbind(col1, col2)))
    
    theta_get = thetaTSMLE(x,
                           lower = -30,
                           upper = 30,
                           p1 = results[[pair[1]]]$p,
                           q1 = results[[pair[1]]]$q,
                           alpha1 = results[[pair[1]]]$alpha,
                           beta1 = results[[pair[1]]]$beta,
                           p2 = results[[pair[2]]]$p,
                           q2 = results[[pair[2]]]$q,
                           alpha2 = results[[pair[2]]]$alpha,
                           beta2 = results[[pair[2]]]$beta
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
    # theta_est = frank_theta_to_spearman(theta_get)
    if(theta_test[1] <=0){
      theta_est = 0
    }
    return(c(pair[1], pair[2], theta_est, theta_test))
  }
  num_cols = ncol(data)
  # Create all combinations of column pairs
  all_pairs <- combn(1:num_cols, 2, simplify = FALSE)
  # Process each chunk
    
    # Set up cluster for this chunk
  cl <- makeCluster(ncores, type = "SOCK")
  registerDoSNOW(cl)
    
    # Progress bar for this chunk
  chunk_pb <- progress_bar$new(
      format = "Chunk :current/:total [:bar] :percent eta: :eta",
      total = length(chunk),
      clear = FALSE,
      width = 60
  )
    
  progress <- function(n) chunk_pb$tick()
  opts <- list(progress = progress)
    
    # Process the chunk
  results_theta <- foreach(pair = all_pairs, 
                             .combine = rbind, 
                             .packages = c('data.table',"copula"), 
                             .options.snow = opts, 
                             .errorhandling = 'pass',
                             .inorder = T,
                             .multicombine = TRUE) %dopar% {
                               compute_theta(pair)
                             }
    
    # Stop the cluster
    # Save results for this chunk
    saveRDS(results_theta, wd_estimate)
    stopCluster(cl)
    
    message("All estimations processed and saved.\n")
    
  
}
