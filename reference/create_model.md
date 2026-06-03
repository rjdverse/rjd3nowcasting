# Create a Dynamic Factor Model

Creates a Dynamic Factor Model (DFM) by specifying the number of factors
and the number of lags in the Vector Auto-Regressive (VAR) process. The
function also allows the user to choose between different types of links
between each series and the latent factors, in addition to defining the
factor loading structure.

## Usage

``` r
create_model(
  nfactors,
  nlags,
  factors_type,
  factors_loading,
  var_init = c("Unconditional", "Zero"),
  var_coefficients = NULL,
  var_errors_variance = NULL,
  measurement_coefficients = NULL,
  measurement_errors_variance = NULL
)
```

## Arguments

- nfactors:

  Integer. Number of factors.

- nlags:

  Integer. Number of lags in the VAR process.

- factors_type:

  Character vector. Defines, in the order of the input series, how each
  (transformed) series relates to the factors. Available options are:

  - "M": Variables expressed as monthly growth rates, linked directly to
    the latent factors.

  - "Q": Quarterly or monthly variables linked to a weighted average of
    the factors, representing quarterly growth rates.

  - "YoY": Variables linked to the cumulative sum of the last 12 monthly
    factors, corresponding to year-on-year growth rates. This option is
    appropriate for variables expressed in year-on-year terms, or for
    series that are closely related to such evolution.

- factors_loading:

  Boolean matrix defining the factor loading structure, with dimension
  *number of series* × *number of factors*. Each entry indicates whether
  a given factor loads on a given series.

- var_init:

  Character. Specifies the initialization of the latent factors: either
  zero or assumed to follow a normal distribution with mean zero and
  variance equal to the unconditional variance of the VAR (default).

- var_coefficients:

  Matrix. The VAR coefficients. If `NULL` (default), they will be
  estimated from scratch. Otherwise, user-defined values can be provided
  as starting point for the estimation, typically obtained from a
  previous estimation of the model.

- var_errors_variance:

  Matrix. The variance-covariance matrix of the VAR errors. If `NULL`
  (default), it will be estimated from scratch. Otherwise, user-defined
  values can be provided as starting point for the estimation, typically
  obtained from a previous estimation of the model.

- measurement_coefficients:

  Matrix. The measurement equation coefficients. If `NULL` (default),
  they will be estimated from scratch. Otherwise, user-defined values
  can be provided as starting point for the estimation, typically
  obtained from a previous estimation of the model.

- measurement_errors_variance:

  Numeric vector. The variance of the idiosyncratic measurement errors.
  If `NULL` (default), they will be estimated from scratch. Otherwise,
  user-defined values can be provided as starting point for the
  estimation, typically obtained from a previous estimation of the
  model.

## Value

An object of class `"JD3_DFMMODEL"` is returned. The following are
returned invisibly as a list:

- `var_coefficients` `[[1]]` initial values of the VAR coefficients;

- `var_errors_variance` `[[2]]` initial values of the
  variance-covariance matrix of the VAR errors;

- `measurement_coefficients` `[[3]]` initial values of the measurement
  equation coefficients;

- `measurement_errors_variance` `[[4]]` initial values of the variance
  of the idiosyncratic measurement errors;

- `initialization_type` `[[5]]` the value of the `var_init` argument;

- `factors_type` `[[6]]` the value of the `factors_type` argument.

## See also

[`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
to estimate the generated model using principal components analysis,

[`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md)
to estimate the generated model using the Expectations-Maximization
algorithm,

[`estimate_ml()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_ml.md)
to estimate the generated model using maximum likelihood.

For more information, see the vignette:

[`utils::browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html),
e.g. `browseVignettes(package = "rjd3nowcasting")`

## Examples

``` r
if (FALSE) { # rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version
# From scratch
dfm1 <- create_model(nfactors = 2,
                     nlags = 2,
                     factors_type = c("M", "M", "YoY", "M", "Q"),
                     factors_loading = matrix(data = TRUE, 5, 2),
                     var_init = "Unconditional")

# From a previously estimated model
set.seed(100)
data <- ts(matrix(rnorm(500), 100, 5),
           frequency = 12,
           start = c(2010, 1))
data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
est1 <- estimate_em(dfm1, data)
dfm2 <- create_model(nfactors = 2,
                     nlags = 2,
                     factors_type = c("M", "M", "YoY", "M", "Q"),
                     factors_loading = matrix(data = TRUE, 5, 2),
                     var_init = "Unconditional",
                     var_coefficients = est1$dfm$var_coefficients,
                     var_errors_variance = est1$dfm$var_errors_variance,
                     measurement_coefficients = est1$dfm$measurement_coefficients,
                     measurement_errors_variance = est1$dfm$measurement_errors_variance)
}
```
