# skewtDist

**skewtDist** is an R package for two families of skew-t distributions that have different tail behavior for left and right tails, namely the family of the Asymmetric Student t-distributions (AST) distributions introduced by Zhu and Galbraith (2010) (ZG), and the the family of generalized asymmetric t-distributions (GAT) introduced by Baker (2018) (BAK). The importance of these two families is that they go beyond the symmetric tail behaviors of the skew-t distributions derived from skew-normal distributions, as described in Azzalini and Capitanio (2014), and hence can provide better fits for certain data arising in applications, especially for asset returns.

The **skewtDist** package aims to provide a set of tools for modelling non-normally distributed data with skew-t distributions. The package includes, for both AST and GAT families, the following functionality: (a) computation of probability density functions, cumulative distribution functions, quantiles, moments (mean, variance, skewness and excess kurtosis) (b) generation of random variables (c) information matrix computation, and (d) maximum likelihood estimation skew-t distribution parameters. 

## Installation

The package is currently in the process of submitting to CRAN.
<!---You can install the released version of st from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("skewtDist")
```--->

``` r
devtools::install__github("dan9401/skewtDist")
```

## Example

``` r
pars <- c(0.12, 0.6, 0.6, 3, 5)
data <- rast(1000, pars = pars)
hist(data, breaks = 50, probability = TRUE)
x <- seq(-3, 3, 0.01)
y <- dast(x, pars = pars)
lines(x, y, col = 4)
```

An example fitting AST distribution using MLE:

``` r
data(retSW)
plxs <- retSW['2006-01/', 'PLXS']
fit <- astMLE(plxs)
summary(fit)
plot(fit, 2)
```

For more examples, please refer to the [Vignette](vignettes/VignetteSkewtDist.pdf).