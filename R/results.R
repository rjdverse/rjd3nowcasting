#' @include utils.R
#' @importFrom stats ts cycle time window
NULL

#' @title Get Dynamic Factor Model Results
#'
#' @description
#' Provides access to estimation results, including pre-processing (e.g.,
#' standardization input for dynamic workflows), parameter and factor estimates,
#' residuals, and likelihood.
#'
#' @param dfm_estimates An object of class `"JD3_DFMESTIMATES"`, typically generated using the `estimate_ml()`, `estimate_em()`, or `estimate_pca()` function.
#'
#' @return An object of class `"JD3_DFMRESULTS"` is returned. The following are
#' returned invisibly as a list:
#' * `preprocessing` `[[1]]` the original and transformed data, along with the sample mean and standard deviation used for standardization;
#' * `parameters` `[[2]]` the estimated parameters, the variance of the measurement errors, and the variance-covariance matrix of the VAR errors;
#' * `factors` `[[3]]` the estimated factors and their standard deviations;
#' * `residuals` `[[4]]` the residuals and the standardized residuals for diagnostic checks;
#' * `likelihood` `[[5]]` the estimated log-likelihood and a Boolean indicating whether the estimation has converged.
#'
#' @export
#'
#' @seealso `get_forecasts()` to obtain forecast results.
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5),
#'            frequency = 12,
#'            start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#'
#' dfm <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' est_em <- estimate_em(dfm, data)
#'
#' rslt_em <- get_results(est_em)
#'
get_results <- function(dfm_estimates) {
    data <- dfm_estimates$data
    dfm <- dfm_estimates$dfm
    is_standardized <- dfm_estimates$is_standardized
    input_standardization <- dfm_estimates$input_standardization
    dfm_rslts <- .dfmProcess(dfm, data, is_standardized, input_standardization)

    return(structure(
        list(
            preprocessing = .get_preprocessing(dfm_rslts),
            parameters = .get_parameters(dfm, dfm_rslts$series_names),
            factors = .get_factors(dfm_rslts),
            residuals = .get_residuals(dfm_rslts),
            likelihood = .get_likelihood(dfm_estimates)
        ),
        class = "JD3_DFMRESULTS"
    ))
}


#' @title Get Dynamic Factor Model Forecasts
#'
#' @description
#' Provides access to forecasts and their associated standard deviations for
#' both the original and transformed series. Optional argument allows to include
#' the missing values estimate in the output.
#'
#' @param dfm_estimates An object of class `"JD3_DFMESTIMATES"`, typically generated using the `estimate_ml()`, `estimate_em()`, or `estimate_pca()` function.
#' @param n_fcst Integer. Number of forecast periods to consider. The default is `3`.
#' @param estim_missing Boolean. Indicates whether missing values should be estimated prior to the start of the forecasting period. The default is `FALSE`.
#' @param mask_q_m Boolean. Indicates whether estimates for the first two months of each quarter should be masked when the factor type is `Q`. The default is `FALSE`.
#'
#' @return An object of class `"JD3_DFMFORECASTS"` is returned. The following are
#'   returned invisibly as a list:
#' * `transformed_forecasts` `[[1]]` the transformed series together with their forecasts;
#' * `transformed_forecasts_stdev` `[[2]]` standard deviations of the transformed series and their forecasts;
#' * `forecasts` `[[3]]` the original series together with their forecasts;
#' * `forecasts_stdev` `[[4]]` standard deviations of the original series and their forecasts;
#' * `forecasts_only` `[[5]]` the forecasts of the original series;
#' * `forecasts_only_stdev` `[[6]]` standard deviations of the forecasts of the original series.
#'
#' @export
#'
#'
#' @seealso `get_results()` to obtain estimation results.
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5),
#'            frequency = 12,
#'            start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#'
#' dfm <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' est_em <- estimate_em(dfm, data)
#'
#' fcsts_em <- get_forecasts(est_em, n_fcst = 2)
#'
get_forecasts <- function(dfm_estimates,
                          n_fcst = 3,
                          estim_missing = FALSE,
                          mask_q_m = FALSE) {

    n_fcst <- max(1, n_fcst)

    data <- dfm_estimates$data
    dfm <- dfm_estimates$dfm
    is_standardized <- dfm_estimates$is_standardized
    input_standardization <- dfm_estimates$input_standardization

    dfm_rslts <- .dfmProcessExt(
        dfm, data, is_standardized, input_standardization, n_fcst
    )

    # Transformed series
    freq <- dfm_rslts$freq
    start <- dfm_rslts$start

    fcsts_t <- stats::ts(
        if (estim_missing) dfm_rslts$forecasts_T_M else dfm_rslts$forecasts_T,
        frequency = freq, start = start
    )
    fcsts_t_stderr <- stats::ts(
        dfm_rslts$forecasts_T_stderr,
        frequency = freq, start = start
    )

    # Original series
    fcsts <- stats::ts(
        if (estim_missing) dfm_rslts$forecasts_M else dfm_rslts$forecasts,
        frequency = freq, start = start
    )
    fcsts_stderr <- stats::ts(
        dfm_rslts$forecasts_stderr,
        frequency = freq, start = start
    )

    colnames(fcsts_t) <- colnames(fcsts_t_stderr) <-
        colnames(fcsts) <- colnames(fcsts_stderr) <-
        dfm_rslts$series_names

    # Mask non-relevant monthly observations for quarterly variables
    if (mask_q_m) {
        f_type <- dfm$factors_type
        period_to_mask <- stats::cycle(fcsts) %% 3 != 0

        for (j in seq_len(ncol(fcsts))) {
            if (f_type[j] == "Q") {
                fcsts[period_to_mask, j] <- NA
            }
        }
    }

    # Restrict output to forecasts only
    nc <- ncol(data)
    nf <- sapply(
        seq_len(nc),
        function(j) {
            cj <- data[, j]
            min(which(!is.na(rev(cj)))) - 1
        }
    )

    nf_max <- max(nf)
    strt <- stats::time(fcsts)[nrow(data) - nf_max + 1]
    strt_yr <- floor(strt)
    strt_mth <- round((strt %% 1) * freq + 1, 0)

    fcsts_only <- stats::window(fcsts, start = c(strt_yr, strt_mth))
    fcsts_only_stderr <- stats::window(fcsts_stderr, start = c(strt_yr, strt_mth))

    ## set observed periods to NA to keep only forecasts
    if (nf_max > 0) {
        for (j in seq_len(nc)) {
            n_na <- nf_max - nf[j]
            if (n_na > 0) {
                fcsts_only[seq_len(n_na), j] <- NA
                fcsts_only_stderr[seq_len(n_na), j] <- NA
            }
        }
    }

    return(structure(
        list(
            transformed_forecasts = fcsts_t,
            transformed_forecasts_stdev = fcsts_t_stderr,
            forecasts = fcsts,
            forecasts_stdev = fcsts_stderr,
            forecasts_only = fcsts_only,
            forecasts_only_stdev = fcsts_only_stderr
        ),
        class = "JD3_DFMFORECASTS"
    ))
}
