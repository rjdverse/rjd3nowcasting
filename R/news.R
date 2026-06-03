#' @include utils.R
#' @importFrom stats frequency start
NULL

#' @title Perform News Analysis for Forecast Updates in Dynamic Factor Models
#'
#' @description
#' Performs a news analysis to quantify the impact of newly released data on
#' forecast revisions between two consecutive dataset updates. The function
#' isolates the contribution of new observations by comparing them to their
#' forecasts obtained from the revised dataset (including past data revisions)
#' and provides a detailed decomposition of the resulting forecast changes.
#'
#' @param dfm_estimates An object of class `"JD3_DFMESTIMATES"`, typically generated using the `estimate_ml()`, `estimate_em()`, or `estimate_pca()` function.
#' @param new_data A `"mts"` object containing the updated dataset.
#' @param target_series Character. Name of the target series of interest. By default, the first series is used.
#' @param n_fcst Integer. Number of forecast periods to consider. The default is `3`.
#'
#' @return An object of class `"JD3_DFMNEWS"` is returned. The following are
#'   returned invisibly as a list:
#' * `target_series` `[[1]]` the value of the `target_series` argument;
#' * `weights_T` `[[2]]` the weights of the news for the transformed series;
#' * `impacts_T` `[[3]]` the impacts of the news (weights x news) for the transformed series;
#' * `forecasts_T` `[[4]]` the old, revised (including past data revisions) and new forecasts of the transformed series;
#' * `weights` `[[5]]` the weights of the news for the original series;
#' * `impacts` `[[6]]` the impacts of the news (weights x news) for the original series;
#' * `forecasts` `[[7]]` the old, revised (including past data revisions) and new forecasts of the original series.
#'
#' @export
#'
#' @seealso `create_model()` to define a new Dynamic Factor model,
#'
#' `estimate_pca()` for estimation using principal components analysis,
#'
#' `estimate_em()` for estimation using EM algorithm,
#'
#' `estimate_ml()` for estimation using maximum likelihood,
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @references
#' Banbura, M., & Modugno, M. (2010). Maximum likelihood estimation of factor
#' models on data sets with arbitrary patterns of missing data.
#'
#' @examplesIf rjd3toolkit::get_java_version() >= rjd3toolkit::minimal_java_version
#' set.seed(100)
#' data_t1 <- ts(matrix(rnorm(500), 100, 5),
#'               frequency = 12,
#'               start = c(2010, 1))
#' data_t1[100, 1] <- data_t1[99:100, 2] <- data_t1[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' data_t2 <- ts(rbind(data_t1, rep(NA, 5)),
#'               frequency = 12,
#'               start = c(2010, 1))
#' data_t2[100, 1] <- data_t2[99, 2] <- data_t2[101, 3] <- data_t2[101, 4] <- 1
#'
#' dfm_model <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' est_em_a <- estimate_em(dfm_model, data_t1)
#'
#' news_a <- get_news(est_em_a, data_t2, target_series = "Series 2", n_fcst = 2)
#'
get_news <- function(dfm_estimates,
                     new_data,
                     target_series = NULL,
                     n_fcst = 3) {

    old_data <- dfm_estimates$data
    series_names <- colnames(old_data)

    if (is.null(target_series)) {
        target_series <- series_names[1]
        target_series_index <- 0

    } else {
        idx <- which(colnames(old_data) == target_series)

        if (length(idx) != 1) {
            warning("Target series not found or not unique. Using the first series instead.")
            target_series <- series_names[1]
            target_series_index <- 0
        } else {
            target_series_index <- idx - 1
        }
    }

    if (n_fcst < 1) n_fcst <- 3
    target_factor_type <- dfm_estimates$dfm$factors_type[target_series_index + 1]
    if (target_factor_type == "Q") {
        n_fcst <- n_fcst * 3
    }

    jdfm <- .r2jd_dfm(dfm_estimates$dfm)
    jold_data <- rjd3toolkit::.r2jd_matrix(old_data)
    jnew_data <- rjd3toolkit::.r2jd_matrix(new_data)
    freq <- stats::frequency(old_data)
    start <- stats::start(old_data)
    jnews <- rjd3toolkit::.jd3_object(
        .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmResultsNews;",
            "computeNews",
            as.integer(target_series_index),
            jdfm$internal,
            jold_data,
            jnew_data,
            as.integer(freq),
            .jarray(as.integer(start)),
            dfm_estimates$is_standardized,
            as.integer(n_fcst)
        ),
        result = TRUE
    )

    rnews_T <- .jd2r_dfmNews(jnews, series_names, transformed = TRUE, target_factor_type)
    rnews <- .jd2r_dfmNews(jnews, series_names, transformed = FALSE, target_factor_type)

    return(structure(
        list(
            target_series = target_series,
            weights_T = rnews_T$weights,
            impacts_T = rnews_T$impacts,
            forecasts_T = rnews_T$fcsts,
            weights = rnews$weights,
            impacts = rnews$impacts,
            forecasts = rnews$fcsts
        ),
        class = "JD3_DFMNEWS"
    ))
}
