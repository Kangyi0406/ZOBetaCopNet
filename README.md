\title{`scCoNet`: Frank Copulas with Zero-One Beta Mixture Margins for Single-Cell Co-expression Covariation Networks}

\section{About}
\code{scCoNet} is an R package for modeling zero-one-inflated proportion data using Frank copula models with mixture margins. It implements a two-stage maximum likelihood estimation (tsMLE) procedure and a two-stage likelihood ratio test (tsLRT) for estimating the copula dependence parameter, which can then be used to build co-expression covariation networks.

\section{Mathematical Model}
\code{scCoNet} uses a mixture margin copula model for pairs of zero-one-inflated beta distributions for gene data $(X_i, X_j)$:

$$
F(x_i, x_j; \gamma_i, \gamma_j, \theta_{ij}) = C\Big( F_i(x_i; \gamma_i),\; F_j(x_j; \gamma_j);\; \theta_{ij} \Big)
$$

where:
- \(F_i\) and \(F_j\) represent the marginal distribution functions with mixture margins.
- \(\theta_{ij}\) is the copula dependence parameter that captures the covariation between the two features.

This model allows the estimation of dependence parameters that serve as the basis for constructing co-expression covariation networks.

\section{Installation}

To install the latest version of \code{scCoNet} from GitHub, use the following command:

\code{
install.packages("devtools")
devtools::install_github("Kangyi0406/scCoNet")
}

\section{Documentation and Examples}

Once installed, you can access the documentation for any function via \code{?} in R. For example:

\code{
library("scCoNet")
?scCoNet
}

### Example Usage:
After loading the data and covariates, run the main function:

\code{
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

\section{Function Parameters}

\describe{
  \item{\code{data}}{A $n \times p$ data frame of gene expression data with rows as samples and columns as genes.}
  \item{\code{covars}}{A $n \times q$ data frame of covariates (e.g., experimental conditions) with rows in the same order as \code{data}.}
  \item{\code{ncores}}{The number of cores for parallelization. Default is 1 (no parallel computing).}
  \item{\code{formula.mu}}{A formula for the mean parameter in the marginal regression model.}
  \item{\code{formula.sigma}}{A formula for the dispersion parameter in the marginal regression model.}
  \item{\code{formula.nu}}{A formula for the zero-one inflation probability parameter.}
  \item{\code{sig_level}}{Significance threshold (after Benjamini–Hochberg adjustment) for determining which associations are retained. Default is 0.05.}
  \item{\code{sig_connect}}{If \code{TRUE}, only statistically significant associations are returned. Default is \code{FALSE}, which returns all estimated connections.}
}

\section{Contact}
For bug reports, issues, or suggestions, please contact the maintainer:
- **Kangyi Zhao** via [email](mailto:kaz78@pitt.edu).
