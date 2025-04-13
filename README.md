# `scCoNet`: Frank Copulas with Zero-One Beta Mixture Margins for Covariation Networks

<!-- badges: start -->
<!-- badges: end -->

ZOBetaCopNet is an R package for modeling zero-one-inflated proportion data using Frank copula models with mixture margins. It implements a two-stage maximum likelihood estimation (tsMLE) procedure and a two-stage likelihood ratio test (tsLRT) for estimating the copula dependence parameter, which can then be used to build covariation networks.

## About

ZOBetaCopNet uses a mixture margin copula model for a pair of zero-inflated microbial relative abundances \((X_i, X_j)\):

$$
F(x_i, x_j; \gamma_i, \gamma_j, \theta_{ij}) = C\Big( F_i(x_i; \gamma_i),\; F_j(x_j; \gamma_j);\; \theta_{ij} \Big)
$$

Here, \(F_i\) and \(F_j\) denote the marginal distribution functions with mixture margins, and the copula dependence parameter \(\theta_{ij}\) captures the covariation between the two features. The estimated dependence parameters serve as a robust basis for constructing and analyzing covariation networks.

## Installation

You can install the latest version of `ZOBetaCopNet` from GitHub with:

```r
install.packages("devtools")
devtools::install_github("Kangyi0406/ZOBetaCopNet")
```

## Documentation and Examples

Help documentation for the ZOBetaCopNet package is available in R. After installing the package from GitHub via devtools and loading it with library(), you can access the documentation for any function using ?. For example:
``` r
?dzib
```



The main function, `copnet` takes two data frames:

1. `abd`: a $n \times p$ data frame of relative abundance, with samples as rows and microbes as columns.
2. `covars`: a $n \times q$ data frame of covariates, with samples as rows in the same order as abd.

The function has one other inputs,  `ncores`: the number of cores for parallelization. The default `ncores` is one core, implying no parallel computing.

## Contact

For bug reports, issues, or suggestions, please contact the maintainer:
 Kangyi Zhao via [email](mailto:kaz78@pitt.edu).
