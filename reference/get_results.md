# Get Dynamic Factor Model Results

Provides access to estimation results, including pre-processing (e.g.,
standardization input for dynamic workflows), parameter and factor
estimates, residuals, and likelihood.

## Usage

``` r
get_results(dfm_estimates)
```

## Arguments

- dfm_estimates:

  An object of class `"JD3_DFMESTIMATES"`, typically generated using the
  [`estimate_ml()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_ml.md),
  [`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md),
  or
  [`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
  function.

## Value

An object of class `"JD3_DFMRESULTS"` is returned. The following are
returned invisibly as a list:

- `preprocessing` `[[1]]` the original and transformed data, along with
  the sample mean and standard deviation used for standardization;

- `parameters` `[[2]]` the estimated parameters, the variance of the
  measurement errors, and the variance-covariance matrix of the VAR
  errors;

- `factors` `[[3]]` the estimated factors and their standard deviations;

- `residuals` `[[4]]` the residuals and the standardized residuals for
  diagnostic checks;

- `likelihood` `[[5]]` the estimated log-likelihood and a Boolean
  indicating whether the estimation has converged.

## See also

[`get_forecasts()`](https://rjdverse.github.io/rjd3nowcasting/reference/get_forecasts.md)
to obtain forecast results.

For more information, see the vignette:

[`utils::browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html),
e.g. `browseVignettes(package = "rjd3nowcasting")`

## Examples

``` r
if (FALSE) { # rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version
set.seed(100)
data <- ts(matrix(rnorm(500), 100, 5),
           frequency = 12,
           start = c(2010, 1))
data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA

dfm <- create_model(
    nfactors = 2,
    nlags = 2,
    factors_type = c("M", "M", "YoY", "M", "Q"),
    factors_loading = matrix(data = TRUE, 5, 2),
    var_init = "Unconditional"
)

est_em <- estimate_em(dfm, data)

rslt_em <- get_results(est_em)
}
```
