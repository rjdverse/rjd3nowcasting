# Estimate a Dynamic Factor Model by Maximum Likelihood

Estimate the model parameters by Maximum Likelihood (ML), with initial
values obtained from the Expectations-Maximization (EM) algorithm by
default. The function includes optional arguments to tune the estimation
process.

## Usage

``` r
estimate_ml(
  dfm,
  data,
  standardized = FALSE,
  input_standardization = NULL,
  pca_init = TRUE,
  em_init = TRUE,
  em_max_iter = 100,
  em_eps = 1e-09,
  max_iter = 1000,
  max_block_iter = 5,
  simpl_model_iter = 15,
  independent_var_shocks = FALSE,
  mixedEstimation = TRUE,
  eps = 1e-09,
  re_estimate = TRUE
)
```

## Arguments

- dfm:

  An object of class `"JD3_DFMMODEL"`, typically generated using the
  [`create_model()`](https://rjdverse.github.io/rjd3nowcasting/reference/create_model.md)
  function.

- data:

  A `"mts"` object containing the (transformed) input data.

- standardized:

  Boolean. Indicates whether the input data are already standardized.
  The default is `FALSE`, in which case standardization is applied as a
  preprocessing step.

- input_standardization:

  A matrix specifying the mean and standard deviation of the series to
  be used during standardization. The default is `NULL`, meaning that
  these quantities are computed from the data. This argument can be set
  using the output of `get_results()$preprocessing$sample_mean_stdev`
  from a previous model estimation. If provided manually, it must be a
  two-column matrix with the means in the first column and the standard
  deviations in the second column. The ordering of the rows must match
  the series in `data`. This argument must be provided if
  `re_estimate = FALSE`, and is ignored if `standardized = TRUE`.

- pca_init:

  Boolean. Indicates whether a Principal Components Analysis (PCA) is
  performed beforehand and used to initialize either the EM algorithm
  (if `em_init = TRUE`) or directly the ML estimation. The default is
  `TRUE`.

- em_init:

  Boolean. Indicates whether the EM algorithm is run prior to ML
  estimation and used to provide initial values. The default is `TRUE`.

- em_max_iter:

  Integer. Specifies the maximum number of iterations for the EM
  algorithm. Ignored if `em_init = FALSE`. The default is `100`.

- em_eps:

  Numeric. The EM algorithm runs until the increase of the percentage
  likelihood falls below the `eps` value (default is `1e-9`), or when
  the maximum number of iterations is reached. Ignored if
  `em_init = FALSE`.

- max_iter:

  Integer. Specifies the maximum number of iterations for the ML
  estimation. The default is `1000`.

- max_block_iter:

  Integer. Specifies the maximum number of iterations per block in the
  optimization process. Model parameters are divided into two blocks:
  one for measurement equation and one for VAR equations. Unlike the EM
  algorithm (one iteration per block), numerical optimization allows
  multiple iterations per block. The default is `5`.

- simpl_model_iter:

  Integer. Specifies the number of iterations allowed for the simplified
  model. The default is `15`.

- independent_var_shocks:

  Boolean. Indicates whether shocks in the VAR block are assumed to be
  independent. The default is `FALSE`.

- mixedEstimation:

  Boolean. Indicates whether to alternate between iterations on the VAR
  block alone and simultaneous iterations on the two blocks. The default
  is `TRUE`.

- eps:

  Numeric. The ML estimation runs until the increase of the percentage
  likelihood falls below the `eps` value (default is `1e-9`), or when
  the maximum number of iterations is reached.

- re_estimate:

  Boolean. Indicates whether the model parameters should be
  re-estimated. The default is `TRUE`. It can be set to `FALSE` to keep
  the model frozen for some period of time, although prolonged use of a
  frozen model is not recommended.

## Value

An object of class `"JD3_DFMESTIMATES"` is returned. The following are
returned invisibly as a list:

- `dfm` `[[1]]` an object of class `"JD3_DFMMODEL"` containing the
  estimated model parameters;

- `data` `[[2]]` the value of the `data` argument;

- `is_standardized` `[[3]]` the value of the `standardized` argument;

- `input_standardization` `[[4]]` the value of the
  `input_standardization` argument;

- `log_likelihood` `[[5]]` the estimated log-likelihood;

- `gradient` `[[6]]` the estimated gradient at the solution;

- `hessian` `[[7]]` the estimated hessian matrix at the solution;

- `has_converged` `[[8]]` a boolean indicating whether the ML estimation
  process has converged.

## See also

[`create_model()`](https://rjdverse.github.io/rjd3nowcasting/reference/create_model.md)
to define a new model,

[`estimate_pca()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_pca.md)
for estimation using principal components analysis,

[`estimate_em()`](https://rjdverse.github.io/rjd3nowcasting/reference/estimate_em.md)
for estimation using EM algorithm,

[`get_results()`](https://rjdverse.github.io/rjd3nowcasting/reference/get_results.md)
to access estimation results,

[`get_forecasts()`](https://rjdverse.github.io/rjd3nowcasting/reference/get_forecasts.md)
to obtain forecasts.

For more information, see the vignette:

[`utils::browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html),
e.g. `browseVignettes(package = "rjd3nowcasting")`

## Examples

``` r
if (FALSE) { # rjd3jars::check_java_version(silent = TRUE)
# input data
set.seed(100)
data <- ts(matrix(rnorm(500), 100, 5),
           frequency = 12,
           start = c(2010, 1))
data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA

# define a new model
dfm <- create_model(
    nfactors = 2,
    nlags = 2,
    factors_type = c("M", "M", "YoY", "M", "Q"),
    factors_loading = matrix(data = TRUE, 5, 2),
    var_init = "Unconditional"
)

# estimate by maximum likelihood
est_ml_a <- estimate_ml(dfm, data)

# results and forecasts
rslts_a <- get_results(est_ml_a)
fcsts_a <- get_forecasts(est_ml_a)

# no re-estimation of a previous model while integrating updated data in the output object
data_new <- data
data_new[99, 2] <- 1
est_ml_b <- estimate_ml(est_ml_a$dfm,
                        data_new,
                        input_standardization = rslts_a$preprocessing$sample_mean_stdev,
                        re_estimate = FALSE)
}
```
