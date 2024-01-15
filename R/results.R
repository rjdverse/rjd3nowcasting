#' Get estimated parameters
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @return a list
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_parameters(rslt_em)
#'
get_parameters <- function(results) {
    # VAR parameters
    pv <- rjd3toolkit::result(results$jestimates, "parameters_var")
    nf <- nrow(pv)
    nl <- ncol(pv) / nf
    pv_cnames <- pv_rnames <- character()
    k <- 1
    for (i in seq_len(nl)) {
        for (j in seq_len(nf)) {
            pv_cnames[k] <- paste0("F", j, "[", -i, "]")
            k <- k + 1
        }
    }
    for (i in seq_len(nf)) {
        pv_rnames[i] <- paste0("F", i)
    }
    colnames(pv) <- pv_cnames
    rownames(pv) <- pv_rnames

    # VAR parameters var-cov
    pvv <- rjd3toolkit::result(results$jestimates, "parameters_var_variance")
    colnames(pvv) <- rownames(pvv) <- pv_rnames

    # Factors parameters
    pf <- rjd3toolkit::result(results$jestimates, "parameters_factors")
    colnames(pf) <- pv_rnames
    rownames(pf) <- results$series_names

    # Factors parameters variance
    pfv <- as.matrix(rjd3toolkit::result(results$jestimates, "parameters_factors_variance"), ncol = 1)
    colnames(pfv) <- "idiosyncratic_variance"
    rownames(pfv) <- results$series_names

    return(list(
        var = pv,
        var_variance = pvv,
        factors = pf,
        factors_variance = pfv
    ))
}

#' Get information about the preprocessing and the input data
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @return a list
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_preprocessing(rslt_em)
#'
get_preprocessing <- function(results) {
    # Sample mean and standard deviation
    smsd <- cbind(
        rjd3toolkit::result(results$jestimates, "sample_mean"),
        rjd3toolkit::result(results$jestimates, "sample_stddev")
    )
    colnames(smsd) <- c("sample_mean", "sample_stddev")
    rownames(smsd) <- results$series_names

    # Input data
    dt <- ts(rjd3toolkit::result(results$jestimates, "input"),
        frequency = results$freq,
        start = results$start
    )
    colnames(dt) <- results$series_names

    # Transformed input data
    dtt <- ts(rjd3toolkit::result(results$jestimates, "input_transformed"),
        frequency = results$freq,
        start = results$start
    )
    colnames(dtt) <- results$series_names

    return(list(
        sample_mean_stdev = smsd,
        data = dt,
        transformed_data = dtt
    ))
}

#' Get estimates of the factors
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @return a list
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_factors(rslt_em)
#'
get_factors <- function(results) {
    # Factors
    f <- ts(rjd3toolkit::result(results$jestimates, "factors"),
        frequency = results$freq,
        start = results$start
    )
    f_cnames <- character()
    for (i in seq_len(ncol(f))) {
        f_cnames[i] <- paste0("F", i)
    }
    colnames(f) <- f_cnames

    # Factors stdev
    fs <- ts(rjd3toolkit::result(results$jestimates, "factors_stderr"),
        frequency = results$freq,
        start = results$start
    )
    colnames(fs) <- f_cnames

    return(list(
        factors = f,
        factors_stdev = fs
    ))
}

#' Get residuals of standardized residuals
#'
#' Residuals are the one-step forecast errors, while the standardized residuals
#' are the one-step forecast errors after standardization (i.e. divided by
#' their standard deviation). The basic diagnostics for normality,
#' heteroscedasticity and auto-correlation should be performed on the
#' standardized residuals.
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @return a list
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_residuals(rslt_em)
#'
get_residuals <- function(results) {
    r <- ts(rjd3toolkit::result(results$jestimates, "residuals"),
        frequency = results$freq,
        start = results$start
    )
    rs <- ts(rjd3toolkit::result(results$jestimates, "residuals_standardized"),
        frequency = results$freq,
        start = results$start
    )
    colnames(r) <- colnames(rs) <- results$series_names

    return(list(
        residuals = r,
        residuals_standardized = rs
    ))
}

#' Get log-likelihood
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @return a numeric value
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_loglikelihood(rslt_em)
#'
get_loglikelihood <- function(results) {
    return(rjd3toolkit::result(results$jestimates, "likelihood_ll"))
}


#' Get forecasts
#'
#' @param results an object of class 'JD3_SsfModelEstimation'
#' @param nf Integer. Number of forecast periods. Must be inferior or equals to
#'   the value of the parameter 'n_out' in the estimation functions.
#' @param forecasts_only Boolean. Whether it returns the whole series or the
#'   forecasts only.
#' @return a list
#' @export
#'
#' @examples
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#' dfm_model <- model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#' rslt_em <- estimate_em(dfm_model, data)
#' get_forecasts(rslt_em, nf = 3, forecasts_only = TRUE)
#'
get_forecasts <- function(results, nf = 12, forecasts_only = TRUE) {
    # Forecasts on original series
    extrct_f <- paste0("forecasts(", nf, ")")
    f <- ts(rjd3toolkit::result(results$jestimates, extrct_f),
        frequency = results$freq,
        start = results$start
    )
    colnames(f) <- results$series_names

    # Standard deviation of the forecasts on original series
    extrct_fs <- paste0("forecasts_stderr(", nf, ")")
    fs <- ts(rjd3toolkit::result(results$jestimates, extrct_fs),
        frequency = results$freq,
        start = results$start
    )
    colnames(fs) <- results$series_names

    # Forecasts on transformed series
    extrct_ft <- paste0("forecasts_transformed(", nf, ")")
    ft <- ts(rjd3toolkit::result(results$jestimates, extrct_ft),
        frequency = results$freq,
        start = results$start
    )
    colnames(ft) <- results$series_names

    # Standard deviation of the forecasts on transformed series
    extrct_fts <- paste0("forecasts_transformed_stderr(", nf, ")")
    fts <- ts(rjd3toolkit::result(results$jestimates, extrct_fts),
        frequency = results$freq,
        start = results$start
    )
    colnames(fts) <- results$series_names

    # Select forecasts only if required
    if (forecasts_only) {
        input <- ts(rjd3toolkit::result(results$jestimates, "input"),
            frequency = results$freq,
            start = results$start
        )

        nc <- ncol(input)
        nf <- vector(mode = "integer", length = nc)
        for (j in seq_len(nc)) {
            cj <- input[, j]
            nf[j] <- min(which(!is.na(rev(cj)))) - 1
        }

        nf_max <- max(nf)
        strt <- time(f)[nrow(input) - nf_max + 1]
        strt_yr <- floor(strt)
        strt_mth <- round((strt %% 1) * results$freq + 1, 0)

        f_f <- window(f, start = c(strt_yr, strt_mth))
        fs_f <- window(fs, start = c(strt_yr, strt_mth))
        ft_f <- window(ft, start = c(strt_yr, strt_mth))
        fts_f <- window(fts, start = c(strt_yr, strt_mth))

        if (nf_max > 0) {
            for (j in seq_len(nc)) {
                n_na <- nf_max - nf[j]
                if (n_na > 0) f_f[seq_len(n_na), j] <- fs_f[seq_len(n_na), j] <- ft_f[seq_len(n_na), j] <- fts_f[seq_len(n_na), j] <- NA
            }
        }
    } else {
        f_f <- f
        fs_f <- fs
        ft_f <- ft
        fts_f <- fts
    }

    return(list(
        forecasts = f_f,
        forecasts_stdev = fs_f,
        forecasts_transformed = ft_f,
        forecasts_transformed_stdev = fts_f
    ))
}

#' Print function for objects of class 'JD3_SsfModelEstimation'
#'
#' @param x an object of class 'JD3_SsfModelEstimation'
#' @export
#'
print.JD3_SsfModelEstimation <- function(x) {
    preproc <- get_preprocessing(x)
    params <- get_parameters(x)

    t1 <- round(cbind(preproc$sample_mean_stdev, params$factors, params$factors_variance), 5)
    factors_name <- paste0("Coeff. of normalized factor F", seq_len(ncol(params$factors)))
    colnames(t1) <- c("Sample mean", "Sample Stdev", factors_name, "Idiosyncratic variance")
    t2 <- round(params$var, 5)
    t3 <- round(params$var_variance, 5)

    cat("Measurement:\n")
    print(t1)
    cat("\n")

    cat("State:\n")
    cat("VAR coefficients:\n")
    print(t2)
    cat("\n")
    cat("Innovative variance:\n")
    print(t3)
}

#' Summary function for objects of class 'JD3_SsfModelEstimation'
#'
#' @param x an object of class 'JD3_SsfModelEstimation'
#' @export
#'
summary.JD3_SsfModelEstimation <- function(x) {
    fcst <- get_forecasts(x, nf = 1, forecasts_only = TRUE)$forecasts
    t1 <- ts(head(fcst, -1), start = start(fcst), frequency = frequency(fcst))

    cat("Nowcasted values (only):\n")
    print(t1)
}

#' Plot function for objects of class 'JD3_SsfModelEstimation'
#'
#' @param x an object of class 'JD3_SsfModelEstimation'
#' @param series_name Character. Name of the series to plot. By default, the
#'   first series will be plotted.
#' @export
#'
plot.JD3_SsfModelEstimation <- function(x, series_name = NULL) {
    fcst_all <- get_forecasts(x, nf = 12, forecasts_only = FALSE)
    fcst <- fcst_all$forecasts
    fcst_stdev <- fcst_all$forecasts_stdev
    fcst_only <- get_forecasts(x, nf = 12, forecasts_only = TRUE)$forecasts

    if (is.null(series_name)) {
        series_name <- colnames(fcst)[1]
    }

    if (series_name %in% colnames(fcst)) {
        s <- fcst[, series_name]
        s_lb <- s - 1.96 * fcst_stdev[, series_name]
        s_ub <- s + 1.96 * fcst_stdev[, series_name]
        sf <- fcst_only[, series_name]
    } else {
        stop("series name not found!")
    }

    ts.plot(s_lb, s_ub, s, sf,
        gpars = list(main = series_name, sub = "Forecasts with 95% CI", xlab = "", ylab = "", lty = c(3, 3, 1, 1), xaxt = "n", type = "o", pch = 20, cex = 0.8, las = 2, col = c("orange", "orange", "black", "red"))
    )
    axis(1, at = seq(start(s)[1], end(s)[1], by = 1), las = 2)
}
