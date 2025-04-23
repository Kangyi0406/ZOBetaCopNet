# `scCoNet`: Frank Copulas with Zero-One Beta Mixture Margins for Single-Cell Co-expression Covariation Networks

<!-- badges: start -->
<!-- badges: end -->

## About

`scCoNet` is an R package for modeling zero-one-inflated proportion data using Frank copula models with mixture margins. It implements a two-stage maximum likelihood estimation (tsMLE) procedure and a two-stage likelihood ratio test (tsLRT) for estimating the copula dependence parameter, which can then be used to build co-expression covariation networks.


`scCoNet` uses a mixture margin copula model for pairs of zero-one-inflated beta distributions for gene data $(X_i, X_j)$:


$$
F(x_i, x_j; \gamma_i, \gamma_j, \theta_{ij}) = C\Big( F_i(x_i; \gamma_i),\; F_j(x_j; \gamma_j);\; \theta_{ij} \Big)
$$

where:

- $F_i$ and $F_j$ represent the marginal distribution functions with mixture margins.
- $\theta_{ij}$ is the copula dependence parameter that captures the covariation between the two features.

This model allows the estimation of dependence parameters that serve as the basis for constructing co-expression covariation networks.


## Installation

To install the latest version of `scCoNet` from GitHub, use the following command:


```r
install.packages("devtools")
devtools::install_github("Kangyi0406/scCoNet")
```

## Documentation and Examples

Once installed, you can access the documentation for any function via `?` in R. For example:

``` r
library("scCoNet")
?scCoNet

```

```{r}
# Load your data (gene values for each cell or feature)
data(gene_data)
data <- gene_data$gen_data

# Load covariates (e.g., experimental conditions or other metadata)
cov <- gene_data$cov_data
cov$name_data <- factor(cov$name_data)
cov$year <- as.numeric(cov$year)
cov <- cov[, 2:3]  # Select the necessary columns

# Run the scCoNet function
result <- scCoNet(
  data,
  covars = data.frame(cov),
  ncores = 1,
  formula.mu = y ~ name_data + year,
  formula.sigma = y ~ name_data + year,
  formula.nu = y ~ name_data + year,
  sig_level = 0.05,
  sig_connect = TRUE
)
```


### Function Parameters

1. `data`: A $n \times p$ data frame of gene expression data with rows as samples and columns as genes.
2. `covars`: A $n \times q$ data frame of covariates (e.g., experimental conditions) with rows in the same order as `data`
3. `ncores`: The number of cores for parallelization. Default is 1 (no parallel computing).
4. `formula.mu`: A formula for the mean parameter in the marginal regression model.
5. `formula.sigma`A formula for the dispersion parameter in the marginal regression model.
6. `formula.nu` A formula for the zero-one inflation probability parameter.
7. `formula.tau` A formula for the zero-one inflation probability parameter.
8. `sig_level` Significance threshold (after Benjamini–Hochberg adjustment) used to determine which associations are retained. Default is 0.05.
9. `sig_connect` If `TRUE`, only statistically significant associations are returned; otherwise, all estimated connections are included. Default is `FALSE`.


## Contact

For bug reports, issues, or suggestions, please contact the maintainer:
 Kangyi Zhao via [email](mailto:kaz78@pitt.edu).
