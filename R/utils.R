#' @importFrom rJava .jpackage .jcall .jnull .jarray .jevalArray .jcast .jcastToArray .jinstanceof is.jnull .jnew .jclass
#' @importFrom stats frequency start ts
NULL

#' @title Datasets including some French macro-economic variables
#'
#' @description
#' The datasets 'data0' and 'data1' acts as successive releases of
#' macro-economic time series. They contain data on monthly industrial
#' production index (PVI), turnover (TURN), quarterly GDP, as well as
#' business survey data (BS) and other survey data (PMI) for both France
#' and the Eurozone. Those datasets are used to illustrate how one of these
#' variable can be nowcasted using the others using a Dynamic Factor model.
#'
#' @name macroIndicators
#' @keywords datasets
"data0"

#' @rdname macroIndicators
"data1"

.r2jd_dfm <- function(dfm_list) {

    m_coef <- dfm_list$measurement_coefficients
    m_var <- dfm_list$measurement_errors_variance
    var_coef <- dfm_list$var_coefficients
    var_var <- dfm_list$var_errors_variance
    factors_type <- dfm_list$factors_type
    var_init <- dfm_list$initialization_type

    nf <- ncol(m_coef)
    nl <- ncol(var_coef) / nf
    loadings <- ifelse(is.nan(m_coef), FALSE, TRUE)
    jloadings <- rjd3toolkit::.r2jd_matrix(loadings)
    jvar_coef <- rjd3toolkit::.r2jd_matrix(var_coef)
    jvar_var <- rjd3toolkit::.r2jd_matrix(var_var)
    jm_coef <- rjd3toolkit::.r2jd_matrix(m_coef)

    jmodel <- .jcall(
        "jdplus/dfm/base/r/DynamicFactorModels",
        "Ljdplus/dfm/base/core/DynamicFactorModel;",
        "model",
        as.integer(nf),
        as.integer(nl),
        .jarray(as.character(factors_type)),
        jloadings,
        .jnew("java/lang/String", as.character(var_init)),
        jvar_coef,
        jvar_var,
        jm_coef,
        .jarray(as.numeric(m_var))
    )

    rjd3toolkit::.jd3_object(jmodel, "JD3_DFMMODEL", result = TRUE)
}

.jd2r_dfm <- function(jmodel) {

    jmodel <- rjd3toolkit::.jd3_object(jmodel, result = TRUE)

    vc <- rjd3toolkit::result(jmodel, "var_coefficients")
    vv <- rjd3toolkit::result(jmodel, "var_errors_variance")
    mc <- rjd3toolkit::result(jmodel, "measurement_coefficients")
    mv <- rjd3toolkit::result(jmodel, "measurement_errors_variance")
    init <- rjd3toolkit::result(jmodel, "initialization_type")
    ftypes <- rjd3toolkit::result(jmodel, "factors_type")

    ftypes_char <- sapply(
        as.character(ftypes),
        switch,
        "1" = "M",
        "5" = "Q",
        "12" = "YoY",
        USE.NAMES = FALSE
    )

    return(structure(
        list(
            var_coefficients = vc,
            var_errors_variance = vv,
            measurement_coefficients = mc,
            measurement_errors_variance = mv,
            initialization_type = init,
            factors_type = ftypes_char
        ),
        class = "JD3_DFMMODEL"
    ))
}

.dfmProcess <- function(dfm,
                        data,
                        standardized = FALSE,
                        input_standardization = NULL) {

    freq <- stats::frequency(data)
    start <- stats::start(data)
    jdata <- rjd3toolkit::.r2jd_matrix(data)
    if (is.null(input_standardization)) {
        standardization_mean <- standardization_stdev <- .jnull(class = "[D")
    } else {
        standardization_mean <- input_standardization[, 1]
        standardization_stdev <- input_standardization[, 2]
    }

    jrslts <- rjd3toolkit::.jd3_object(
        .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmResults;",
            "process",
            .r2jd_dfm(dfm)$internal,
            jdata,
            as.integer(freq),
            .jarray(as.integer(start)),
            standardized,
            standardization_mean,
            standardization_stdev,
            as.integer(0)
        ),
        result = TRUE
    )

    return(
        list(
            series_names = colnames(data),
            start = start,
            freq = freq,
            original_data = rjd3toolkit::result(jrslts, "input"),
            sample_mean = rjd3toolkit::result(jrslts, "sample_mean"),
            sample_stddev = rjd3toolkit::result(jrslts, "sample_stddev"),
            transformed_data = rjd3toolkit::result(jrslts, "input_transformed"),
            factors = rjd3toolkit::result(jrslts, "factors"),
            factors_stderr = rjd3toolkit::result(jrslts, "factors_stderr"),
            factors = rjd3toolkit::result(jrslts, "factors"),
            residuals = rjd3toolkit::result(jrslts, "residuals"),
            standard_residuals = rjd3toolkit::result(jrslts, "residuals_standardized")
        )
    )
}

# n_out: Number of out-of-sample periods used during processing. Must be greater
# than the number of forecast periods requested in the forecasting function.
.dfmProcessExt <- function(dfm,
                           data,
                           standardized = FALSE,
                           input_standardization = NULL,
                           n_out = 12) {

    freq <- stats::frequency(data)
    start <- stats::start(data)
    jdata <- rjd3toolkit::.r2jd_matrix(data)
    if (is.null(input_standardization)) {
        standardization_mean <- standardization_stdev <- .jnull(class = "[D")
    } else {
        standardization_mean <- input_standardization[, 1]
        standardization_stdev <- input_standardization[, 2]
    }
    jrslts <- rjd3toolkit::.jd3_object(
        .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmResults;",
            "process",
            .r2jd_dfm(dfm)$internal,
            jdata,
            as.integer(freq),
            .jarray(as.integer(start)),
            standardized,
            standardization_mean,
            standardization_stdev,
            as.integer(n_out)
        ),
        result = TRUE
    )

    return(
        list(
            series_names = colnames(data),
            start = start,
            freq = freq,
            original_data = rjd3toolkit::result(jrslts, "input"),
            sample_mean = rjd3toolkit::result(jrslts, "sample_mean"),
            sample_stddev = rjd3toolkit::result(jrslts, "sample_stddev"),
            transformed_data = rjd3toolkit::result(jrslts, "input_transformed"),
            forecasts_T = rjd3toolkit::result(jrslts, paste0("forecasts_transformed(", n_out, ")")),
            forecasts_T_M = rjd3toolkit::result(jrslts, paste0("forecasts_transformed_miss(", n_out, ")")),
            forecasts_T_stderr = rjd3toolkit::result(jrslts, paste0("forecasts_transformed_stderr(", n_out, ")")),
            forecasts = rjd3toolkit::result(jrslts, paste0("forecasts(", n_out, ")")),
            forecasts_M = rjd3toolkit::result(jrslts, paste0("forecasts_miss(", n_out, ")")),
            forecasts_stderr = rjd3toolkit::result(jrslts, paste0("forecasts_stderr(", n_out, ")"))
        )
    )
}

.get_preprocessing <- function(dfm_rslts) {
    series_names <- dfm_rslts$series_names

    # Original data
    data <- stats::ts(dfm_rslts$original_data,
                      frequency = dfm_rslts$freq,
                      start = dfm_rslts$start)
    colnames(data) <- dfm_rslts$series_names

    # Sample mean and standard deviation
    sample_mean_stddev <- cbind(dfm_rslts$sample_mean, dfm_rslts$sample_stddev)
    colnames(sample_mean_stddev) <- c("sample_mean", "sample_stddev")
    rownames(sample_mean_stddev) <- series_names

    # Transformed data
    data_t <- stats::ts(
        dfm_rslts$transformed_data,
        frequency = dfm_rslts$freq,
        start = dfm_rslts$start
    )
    colnames(data_t) <- series_names

    return(
        list(
            original_data = data,
            sample_mean_stdev = sample_mean_stddev,
            transformed_data = data_t
        )
    )
}

.get_parameters <- function(dfm, series_names) {

    # VAR coefficients
    var_coef <- dfm$var_coefficients
    nfactors <- nrow(var_coef)
    nlags <- ncol(var_coef) / nfactors

    var_coef_cnames <- var_coef_rnames <- character()
    k <- 1
    for (i in 1:nlags) {
        for (j in 1:nfactors) {
            if (i == 1) var_coef_rnames[j] <- paste0("F", j)
            var_coef_cnames[k] <- paste0("F", j, "[", -i, "]")
            k <- k + 1
        }
    }
    colnames(var_coef) <- var_coef_cnames
    rownames(var_coef) <- var_coef_rnames

    # Variance-covariance matrix of the VAR errors
    var_err_variance <- dfm$var_errors_variance
    colnames(var_err_variance) <- rownames(var_err_variance) <- var_coef_rnames

    # Measurement equation coefficients
    measurement_coef <- dfm$measurement_coefficients
    colnames(measurement_coef) <- var_coef_rnames
    rownames(measurement_coef) <- series_names

    # Variance of the idiosyncratic measurement errors
    measurement_err_variance <- as.matrix(dfm$measurement_errors_variance, ncol = 1)
    colnames(measurement_err_variance) <- "idiosyncratic_variance"
    rownames(measurement_err_variance) <- series_names

    return(
        list(
            var_coefficients = var_coef,
            var_errors_variance = var_err_variance,
            measurement_coefficients = measurement_coef,
            measurement_errors_variance = measurement_err_variance
        )
    )
}

.get_factors <- function(dfm_rslts) {

    # Factors
    factors <- stats::ts(dfm_rslts$factors,
                         frequency = dfm_rslts$freq,
                         start = dfm_rslts$start)
    factors_cnames <- character()
    for (i in seq_len(ncol(factors))) {
        factors_cnames[i] <- paste0("F", i)
    }
    colnames(factors) <- factors_cnames

    # Standard deviation of the factors
    factors_stdev <- stats::ts(dfm_rslts$factors_stderr,
                               frequency = dfm_rslts$freq,
                               start = dfm_rslts$start)
    colnames(factors_stdev) <- factors_cnames

    return(list(factors = factors, factors_stdev = factors_stdev))
}


# Residuals correspond to one-step-ahead forecast errors, while standardized
# residuals are obtained by dividing these errors by their standard deviation.
# Diagnostic checks for normality, heteroscedasticity, and autocorrelation
# should be performed on the standardized residuals.
.get_residuals <- function(dfm_rslts) {
    res <- stats::ts(dfm_rslts$residuals,
                     frequency = dfm_rslts$freq,
                     start = dfm_rslts$start)
    res_std <- stats::ts(
        dfm_rslts$standard_residuals,
        frequency = dfm_rslts$freq,
        start = dfm_rslts$start
    )
    colnames(res) <- colnames(res_std) <- dfm_rslts$series_names

    return(list(residuals = res, standardized_residuals = res_std))
}

.get_likelihood <- function(dfm_estimates) {
    ll <- dfm_estimates$log_likelihood
    has_converged <- dfm_estimates$has_converged

    return(list(log_likelihood = ll, has_converged = has_converged))
}


.jd2r_dfmNews <- function(jnews,
                          series_names,
                          transformed = TRUE,
                          target_factor_type = "M") {

    suffix <- if (transformed) "_T" else ""

    fcst_periods <- rjd3toolkit::result(jnews, "forecasts_periods")

    base <- data.frame(
        series = series_names[rjd3toolkit::result(jnews, "series_index") + 1],
        period = rjd3toolkit::result(jnews, "series_period"),
        expected_value = rjd3toolkit::result(jnews, paste0("series_expected_value", suffix)),
        observed_value = rjd3toolkit::result(jnews, paste0("series_observed_value", suffix)),
        news = rjd3toolkit::result(jnews, paste0("series_news", suffix))
    )

    # Weights and impacts
    w <- rjd3toolkit::result(jnews, paste0("series_weights", suffix))
    i <- rjd3toolkit::result(jnews, paste0("series_impacts", suffix))

    colnames(w) <- paste0("weights(", fcst_periods, ")")
    colnames(i) <- paste0("impacts(", fcst_periods, ")")

    weights <- cbind(base, w)

    total_row <- data.frame(
        series = "TOTAL",
        period = NA,
        expected_value = NA,
        observed_value = NA,
        news = NA
    )
    total_row <- cbind(total_row, as.data.frame(t(colSums(i))))

    impacts <- rbind(cbind(base, i), total_row)

    # Forecasts
    fcsts <- rbind(
        rjd3toolkit::result(jnews, paste0("old_forecasts", suffix)),
        rjd3toolkit::result(jnews, paste0("revised_forecasts", suffix)),
        rjd3toolkit::result(jnews, paste0("new_forecasts", suffix))
    )
    colnames(fcsts) <- fcst_periods
    rownames(fcsts) <- c("old_forecasts", "revised_forecasts", "new_forecasts")

    ## adjustment for quarterly variables
    if (target_factor_type == "Q") {

        nc <- ncol(fcsts)
        offset <- ncol(base)
        keep <- logical(nc)

        for (j in seq_len(nc)) {
            label <- fcst_periods[j]
            nk <- nchar(label)

            mth <- as.numeric(substr(label, 1, nk - 5))
            yr <- as.numeric(substr(label, nk - 3, nk))

            if (mth %% 3 == 0) {
                q <- mth / 3
                keep[j] <- TRUE

                colnames(fcsts)[j] <- paste0("Q", q, "-", yr)
                colnames(weights)[j + offset] <- paste0("weights(Q", q, "-", yr, ")")
                colnames(impacts)[j + offset] <- paste0("Impacts(Q", q, "-", yr, ")")
            }
        }

        fcsts <- fcsts[, keep, drop = FALSE]
        weights <- weights[, c(seq_len(offset), offset + which(keep)), drop = FALSE]
        impacts <- impacts[, c(seq_len(offset), offset + which(keep)), drop = FALSE]
    }

    # Output
    list(
        weights = weights,
        impacts = impacts,
        fcsts = fcsts
    )
}
