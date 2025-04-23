\name{scCoNet}
\alias{scCoNet}
\title{Frank Copulas with Zero-One Beta Mixture Margins for Single-Cell Co-expression Covariation Networks}
\description{
  \code{scCoNet} is an R package for modeling zero-one-inflated proportion data using Frank copula models with mixture margins. It implements a two-stage maximum likelihood estimation (tsMLE) procedure and a two-stage likelihood ratio test (tsLRT) for estimating the copula dependence parameter, which can then be used to build co-expression covariation networks.
}
\usage{
  scCoNet(data, covars, ncores = 1, formula.mu, formula.sigma, formula.nu, sig_level = 0.05, sig_connect = FALSE)
}
\arguments{
  \item{data}{A data frame of gene expression data (size: $n \times p$), where rows represent samples (cells) and columns represent genes.}
  \item{covars}{A data frame of covariates (size: $n \times q$), where rows represent samples, and columns represent the covariates. The order of rows must match those of \code{data}.}
  \item{ncores}{Number of cores for parallelization. Default is 1, meaning no parallel computing.}
  \item{formula.mu}{A formula for the mean parameter in the marginal regression model.}
  \item{formula.sigma}{A formula for the dispersion parameter in the marginal regression model.}
  \item{formula.nu}{A formula for the zero-one inflation probability parameter.}
  \item{sig_level}{The significance threshold (after Benjamini–Hochberg adjustment) for determining which associations are retained. Default is 0.05.}
  \item{sig_connect}{If \code{TRUE}, only statistically significant associations are returned. If \code{FALSE}, all estimated connections are included. Default is \code{FALSE}.}
}
\value{
  A list containing the estimated dependence parameters (copula) and other relevant information such as significance levels for the associations.
}
\details{
  \code{scCoNet} uses a mixture margin copula model for pairs of zero-one-inflated beta distributions for gene expression data. Specifically, it models the relationship between two genes \(X_i\) and \(X_j\) using the Frank copula, defined as:

  $$
  F(x_i, x_j; \gamma_i, \gamma_j, \theta_{ij}) = C\left( F_i(x_i; \gamma_i), F_j(x_j; \gamma_j); \theta_{ij} \right)
  $$

  where \(F_i\) and \(F_j\) are the marginal distribution functions for the respective genes, and \(\theta_{ij}\) is the copula dependence parameter that captures the covariation between the two features. This parameter serves as the foundation for constructing co-expression covariation networks.
}
\examples{
  # Load your data (gene values for each cell or feature)
  data <- data.frame(readRDS("path/to/data.rds"))

  # Load covariates (e.g., experimental conditions or other metadata)
  cov <- as.data.frame(readRDS("path/to/covariates.rds"))
  cov$name_data <- factor(cov$name_data)
  cov$year <- as.numeric(cov$year)
  cov <- cov[, 2:3]  # Select the necessary columns

  # Run the scCoNet function
  result <- scCoNet(
    data,
    covars = data.frame(cov),
    ncores = 2,
    formula.mu = y ~ name_data + year,
    formula.sigma = y ~ name_data + year,
    formula.nu = y ~ name_data + year,
    sig_level = 0.05,
    sig_connect = TRUE
  )
}
\seealso{
  \code{\link{other_function_name}} for related functions and models used for co-expression network construction.
}
\author{
  Kangyi Zhao <kaz78@pitt.edu>
}
\keywords{copula, zero-one-inflated data, single-cell, co-expression, copula dependence, statistical modeling}
