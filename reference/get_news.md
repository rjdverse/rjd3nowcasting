# Perform News Analysis for Forecast Updates in Dynamic Factor Models

Performs a news analysis to quantify the impact of newly released data
on forecast revisions between two consecutive dataset updates. The
function isolates the contribution of new observations by comparing them
to their forecasts obtained from the revised dataset (including past
data revisions) and provides a detailed decomposition of the resulting
forecast changes.

## Usage

``` r
get_news(dfm_estimates, new_data, target_series = NULL, n_fcst = 3)
```

## Arguments

- dfm_estimates:

  An object of class `"JD3_DFMESTIMATES"`, typically generated using the
  [`estimate_ml()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_ml.md),
  [`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md),
  or
  [`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
  function.

- new_data:

  A `"mts"` object containing the updated dataset.

- target_series:

  Character. Name of the target series of interest. By default, the
  first series is used.

- n_fcst:

  Integer. Number of forecast periods to consider. The default is `3`.

## Value

An object of class `"JD3_DFMNEWS"` is returned. The following are
returned invisibly as a list:

- `target_series` `[[1]]` the value of the `target_series` argument;

- `weights_T` `[[2]]` the weights of the news for the transformed
  series;

- `impacts_T` `[[3]]` the impacts of the news (weights x news) for the
  transformed series;

- `forecasts_T` `[[4]]` the old, revised (including past data revisions)
  and new forecasts of the transformed series;

- `weights` `[[5]]` the weights of the news for the original series;

- `impacts` `[[6]]` the impacts of the news (weights x news) for the
  original series;

- `forecasts` `[[7]]` the old, revised (including past data revisions)
  and new forecasts of the original series.

## References

Banbura, M., & Modugno, M. (2010). Maximum likelihood estimation of
factor models on data sets with arbitrary patterns of missing data.

## See also

[`create_model()`](https://rjdverse.github.io/rjd3nowcasting/reference/create_model.md)
to define a new Dynamic Factor model,

[`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
for estimation using principal components analysis,

[`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md)
for estimation using EM algorithm,

[`estimate_ml()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_ml.md)
for estimation using maximum likelihood,

For more information, see the vignette:

[`utils::browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html),
e.g. `browseVignettes(package = "rjd3nowcasting")`

## Examples

``` r
if (FALSE) { # rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version
set.seed(100)
data_t1 <- ts(matrix(rnorm(500), 100, 5),
              frequency = 12,
              start = c(2010, 1))
data_t1[100, 1] <- data_t1[99:100, 2] <- data_t1[(1:100)[-seq(3, 100, 3)], 5] <- NA
data_t2 <- ts(rbind(data_t1, rep(NA, 5)),
              frequency = 12,
              start = c(2010, 1))
data_t2[100, 1] <- data_t2[99, 2] <- data_t2[101, 3] <- data_t2[101, 4] <- 1

dfm_model <- create_model(
    nfactors = 2,
    nlags = 2,
    factors_type = c("M", "M", "YoY", "M", "Q"),
    factors_loading = matrix(TRUE, 5, 2),
    var_init = "Unconditional"
)

est_em_a <- estimate_em(dfm_model, data_t1)

news_a <- get_news(est_em_a, data_t2, target_series = "Series 2", n_fcst = 2)
}
```
