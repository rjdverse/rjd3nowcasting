# Get Dynamic Factor Model Forecasts

Provides access to forecasts and their associated standard deviations
for both the original and transformed series. Optional argument allows
to include the missing values estimate in the output.

## Usage

``` r
get_forecasts(
  dfm_estimates,
  n_fcst = 3,
  estim_missing = FALSE,
  mask_q_m = FALSE
)
```

## Arguments

- dfm_estimates:

  An object of class `"JD3_DFMESTIMATES"`, typically generated using the
  [`estimate_ml()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_ml.md),
  [`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md),
  or
  [`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
  function.

- n_fcst:

  Integer. Number of forecast periods to consider. The default is `3`.

- estim_missing:

  Boolean. Indicates whether missing values should be estimated prior to
  the start of the forecasting period. The default is `FALSE`.

- mask_q_m:

  Boolean. Indicates whether estimates for the first two months of each
  quarter should be masked when the factor type is `Q`. The default is
  `FALSE`.

## Value

An object of class `"JD3_DFMFORECASTS"` is returned. The following are
returned invisibly as a list:

- `transformed_forecasts` `[[1]]` the transformed series together with
  their forecasts;

- `transformed_forecasts_stdev` `[[2]]` standard deviations of the
  transformed series and their forecasts;

- `forecasts` `[[3]]` the original series together with their forecasts;

- `forecasts_stdev` `[[4]]` standard deviations of the original series
  and their forecasts;

- `forecasts_only` `[[5]]` the forecasts of the original series;

- `forecasts_only_stdev` `[[6]]` standard deviations of the forecasts of
  the original series.

## See also

[`get_results()`](https://rjdverse.github.io/rjd3nowcasting/reference/get_results.md)
to obtain estimation results.

For more information, see the vignette:

[`utils::browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html),
e.g. `browseVignettes(package = "rjd3nowcasting")`

## Examples

``` r
if (FALSE) { # rjd3jars::check_java_version(silent = TRUE)
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

fcsts_em <- get_forecasts(est_em, n_fcst = 2)
}
```
