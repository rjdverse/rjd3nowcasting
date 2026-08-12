#' @include utils.R
#' @importFrom stats frequency start
NULL

#' @title Estimate a Dynamic Factor Model Using Principal Components Analysis
#'
#' @description
#' Estimate the model parameters using Principal Component Analysis (PCA).
#' While this approach is fast, it is generally not recommended to rely on it
#' exclusively, particularly when the dataset includes quarterly series or
#' variables related to year-on-year growth rates.
#'
#' @param dfm An object of class `"JD3_DFMMODEL"`, typically generated using the `create_model()` function.
#' @param data A `"mts"` object containing the (transformed) input data.
#' @param standardized Boolean. Indicates whether the input data are already standardized.
#' The default is `FALSE`, in which case standardization is applied as a preprocessing step.
#' @param input_standardization A matrix specifying the mean and standard deviation of the series to be used during standardization.
#' The default is `NULL`, meaning that these quantities are computed from the data.
#' This argument can be set using the output of `get_results()$preprocessing$sample_mean_stdev` from a previous model estimation. If provided manually, it must be a two-column matrix with the means in the first column and the standard deviations in the second column. The ordering of the rows must match the series in `data`.
#' This argument must be provided if `re_estimate = FALSE`, and is ignored if `standardized = TRUE`.
#' @param re_estimate Boolean. Indicates whether the model parameters should be re-estimated.
#' The default is `TRUE`. It can be set to `FALSE` to keep the model frozen for some period of time, although prolonged use of a frozen model is not recommended.
#'
#' @return An object of class `"JD3_DFMESTIMATES"` is returned. The following are
#'   returned invisibly as a list:
#' * `dfm` `[[1]]` an object of class `"JD3_DFMMODEL"` containing the estimated model parameters;
#' * `data` `[[2]]` the value of the `data` argument;
#' * `is_standardized` `[[3]]` the value of the `standardized` argument;
#' * `input_standardization` `[[4]]` the value of the `input_standardization` argument;
#' * additional elements not relevant for Principal Components Analysis estimates.
#'
#' @export
#'
#' @seealso `create_model()` to define a new model,
#'
#' `estimate_em()` for estimation using the Expectations-Maximization algorithm,
#'
#' `estimate_ml()` for estimation using maximum likelihood,
#'
#' `get_results()` to access estimation results,
#'
#' `get_forecasts()` to obtain forecasts.
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' # input data
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5),
#'            frequency = 12,
#'            start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#'
#' # define a new model
#' dfm <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "M", "M", "M"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' # estimate using PCA
#' est_pca_a <- estimate_pca(dfm, data)
#'
#' # results and forecasts
#' rslts_a <- get_results(est_pca_a)
#' fcsts_a <- get_forecasts(est_pca_a)
#'
#' # no re-estimation of a previous model while integrating updated data in the output object
#' data_new <- data
#' data_new[99, 2] <- 1
#' est_pca_b <- estimate_pca(est_pca_a$dfm,
#'                           data_new,
#'                           input_standardization = rslts_a$preprocessing$sample_mean_stdev,
#'                           re_estimate = FALSE)
#'
estimate_pca <- function(dfm,
                         data,
                         standardized = FALSE,
                         input_standardization = NULL,
                         re_estimate = TRUE) {

    if (re_estimate) {
        jdfm <- .r2jd_dfm(dfm)
        freq <- stats::frequency(data)
        start <- stats::start(data)
        jdata <- rjd3toolkit::.r2jd_matrix(data)
        if (is.null(input_standardization)) {
            standardization_mean <- standardization_stdev <- .jnull(class = "[D")
        } else {
            standardization_mean <- input_standardization[, 1]
            standardization_stdev <- input_standardization[, 2]
        }

        jest <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmEstimates;",
            "estimate_PCA",
            jdfm$internal,
            jdata,
            as.integer(freq),
            .jarray(as.integer(start)),
            standardized,
            standardization_mean,
            standardization_stdev
        )

        jmodel <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DynamicFactorModel;",
            "getDfm",
            jest
        )

        dfm_list <- .jd2r_dfm(jmodel)
        jest <- rjd3toolkit::.jd3_object(jest, result = TRUE)
        ll <- rjd3toolkit::result(jest, "likelihood_ll")
        gradient <- rjd3toolkit::result(jest, "gradient")
        hessian <- rjd3toolkit::result(jest, "hessian")
        has_converged <- rjd3toolkit::result(jest, "has_converged")

    } else {
        if (!standardized && is.null(input_standardization)) {
            stop("Re-estimation is turned off. Please provide the original mean and standard deviation in 'input_standardization'.")
        }

        dfm_list <- dfm
        ll <- NA
        gradient <- numeric()
        hessian <- matrix()
        has_converged <- "not re-estimated"
    }

    return(structure(
        list(
            dfm = dfm_list,
            data = data,
            is_standardized = standardized,
            input_standardization = input_standardization,
            log_likelihood = ll,
            gradient = gradient,
            hessian = hessian,
            has_converged = has_converged
        ),
        class = "JD3_DFMESTIMATES"
    ))
}

#' @title Estimate a Dynamic Factor Model Using the Expectations-Maximization Algorithm
#'
#' @description
#' Estimate the model parameters using the Expectations-Maximization (EM)
#' algorithm, with PCA-based initialization by default. The function includes
#' optional arguments to tune the estimation process.
#'
#' @param dfm An object of class `"JD3_DFMMODEL"`, typically generated using the `create_model()` function.
#' @param data A `"mts"` object containing the (transformed) input data.
#' @param standardized Boolean. Indicates whether the input data are already standardized.
#' The default is `FALSE`, in which case standardization is applied as a preprocessing step.
#' @param input_standardization A matrix specifying the mean and standard deviation of the series to be used during standardization.
#' The default is `NULL`, meaning that these quantities are computed from the data.
#' This argument can be set using the output of `get_results()$preprocessing$sample_mean_stdev` from a previous model estimation. If provided manually, it must be a two-column matrix with the means in the first column and the standard deviations in the second column. The ordering of the rows must match the series in `data`.
#' This argument must be provided if `re_estimate = FALSE`, and is ignored if `standardized = TRUE`.
#' @param pca_init Boolean. Indicates whether a Principal Components Analysis (PCA) is performed beforehand and used to initialize the EM algorithm. The default is `TRUE`.
#' @param max_iter Integer. Specifies the maximum number of iterations in the EM algorithm. The default is `100`.
#' @param eps Numeric. The EM algorithm runs until the increase of the percentage likelihood falls below the `eps` value (default is `1e-9`), or when the maximum number of iterations is reached.
#' @param re_estimate Boolean. Indicates whether the model parameters should be re-estimated.
#' The default is `TRUE`. It can be set to `FALSE` to keep the model frozen for some period of time, although prolonged use of a frozen model is not recommended.
#'
#' @return An object of class `"JD3_DFMESTIMATES"` is returned. The following are
#'   returned invisibly as a list:
#' * `dfm` `[[1]]` an object of class `"JD3_DFMMODEL"` containing the estimated model parameters;
#' * `data` `[[2]]` the value of the `data` argument;
#' * `is_standardized` `[[3]]` the value of the `standardized` argument;
#' * `input_standardization` `[[4]]` the value of the `input_standardization` argument;
#' * `log_likelihood` `[[5]]` the estimated log-likelihood;
#' * a couple of other elements not relevant for the EM estimates;
#' * `has_converged` `[[8]]` a boolean indicating whether the EM algorithm has converged.
#'
#' @export
#'
#' @seealso `create_model()` to define a new model,
#'
#' `estimate_pca()` for estimation using principal components analysis,
#'
#' `estimate_ml()` for estimation using maximum likelihood,
#'
#' `get_results()` to access estimation results,
#'
#' `get_forecasts()` to obtain forecasts.
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' # input data
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5),
#'            frequency = 12,
#'            start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#'
#' # define a new model
#' dfm <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' # estimate using EM algorithm
#' est_em_a <- estimate_em(dfm, data)
#'
#' # results and forecasts
#' rslts_a <- get_results(est_em_a)
#' fcsts_a <- get_forecasts(est_em_a)
#'
#' # no re-estimation of a previous model while integrating updated data in the output object
#' data_new <- data
#' data_new[99, 2] <- 1
#' est_em_b <- estimate_em(est_em_a$dfm,
#'                         data_new,
#'                         input_standardization = rslts_a$preprocessing$sample_mean_stdev,
#'                         re_estimate = FALSE)
#'
estimate_em <- function(dfm,
                        data,
                        standardized = FALSE,
                        input_standardization = NULL,
                        pca_init = TRUE,
                        max_iter = 100,
                        eps = 1e-9,
                        re_estimate = TRUE) {

    if (re_estimate) {
        jdfm <- .r2jd_dfm(dfm)
        freq <- stats::frequency(data)
        start <- stats::start(data)
        jdata <- rjd3toolkit::.r2jd_matrix(data)
        if (is.null(input_standardization)) {
            standardization_mean <- standardization_stdev <- .jnull(class = "[D")
        } else {
            standardization_mean <- input_standardization[, 1]
            standardization_stdev <- input_standardization[, 2]
        }

        jest <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmEstimates;",
            "estimate_EM",
            jdfm$internal,
            jdata,
            as.integer(freq),
            .jarray(as.integer(start)),
            standardized,
            standardization_mean,
            standardization_stdev,
            pca_init,
            as.integer(max_iter),
            as.numeric(eps)
        )

        jmodel <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DynamicFactorModel;",
            "getDfm",
            jest
        )

        dfm_list <- .jd2r_dfm(jmodel)
        jest <- rjd3toolkit::.jd3_object(jest, result = TRUE)
        ll <- rjd3toolkit::result(jest, "likelihood_ll")
        gradient <- rjd3toolkit::result(jest, "gradient")
        hessian <- rjd3toolkit::result(jest, "hessian")
        has_converged <- rjd3toolkit::result(jest, "has_converged")

    } else {
        if (!standardized && is.null(input_standardization)) {
            stop("Re-estimation is turned off. Please provide the original mean and standard deviation in 'input_standardization'.")
        }

        dfm_list <- dfm
        ll <- NA
        gradient <- numeric()
        hessian <- matrix()
        has_converged <- "not re-estimated"
    }

    return(structure(
        list(
            dfm = dfm_list,
            data = data,
            is_standardized = standardized,
            input_standardization = input_standardization,
            log_likelihood = ll,
            gradient = gradient,
            hessian = hessian,
            has_converged = has_converged
        ),
        class = "JD3_DFMESTIMATES"
    )
    )
}


#' @title Estimate a Dynamic Factor Model by Maximum Likelihood
#'
#' @description
#' Estimate the model parameters by Maximum Likelihood (ML), with initial
#' values obtained from the Expectations-Maximization (EM) algorithm by default.
#' The function includes optional arguments to tune the estimation process.
#'
#' @param dfm An object of class `"JD3_DFMMODEL"`, typically generated using the `create_model()` function.
#' @param data A `"mts"` object containing the (transformed) input data.
#' @param standardized Boolean. Indicates whether the input data are already standardized.
#' The default is `FALSE`, in which case standardization is applied as a preprocessing step.
#' @param input_standardization A matrix specifying the mean and standard deviation of the series to be used during standardization.
#' The default is `NULL`, meaning that these quantities are computed from the data.
#' This argument can be set using the output of `get_results()$preprocessing$sample_mean_stdev` from a previous model estimation. If provided manually, it must be a two-column matrix with the means in the first column and the standard deviations in the second column. The ordering of the rows must match the series in `data`.
#' This argument must be provided if `re_estimate = FALSE`, and is ignored if `standardized = TRUE`.
#' @param pca_init Boolean. Indicates whether a Principal Components Analysis (PCA) is performed beforehand and used to initialize either the EM algorithm (if `em_init = TRUE`) or directly the ML estimation. The default is `TRUE`.
#' @param em_init Boolean. Indicates whether the EM algorithm is run prior to ML estimation and used to provide initial values. The default is `TRUE`.
#' @param em_max_iter Integer. Specifies the maximum number of iterations for the EM algorithm. Ignored if `em_init = FALSE`. The default is `100`.
#' @param em_eps Numeric. The EM algorithm runs until the increase of the percentage likelihood falls below the `eps` value (default is `1e-9`), or when the maximum number of iterations is reached. Ignored if `em_init = FALSE`.
#' @param max_iter Integer. Specifies the maximum number of iterations for the ML estimation. The default is `1000`.
#' @param max_block_iter Integer. Specifies the maximum number of iterations per block in the optimization process.
#' Model parameters are divided into two blocks: one for measurement equation and one for VAR equations. Unlike the EM algorithm (one iteration per block), numerical optimization allows multiple iterations per block. The default is `5`.
#' @param simpl_model_iter Integer. Specifies the number of iterations allowed for the simplified model. The default is `15`.
#' @param independent_var_shocks Boolean. Indicates whether shocks in the VAR block are assumed to be independent. The default is `FALSE`.
#' @param mixedEstimation Boolean. Indicates whether to alternate between iterations on the VAR block alone and simultaneous iterations on the two blocks. The default is `TRUE`.
#' @param eps Numeric. The ML estimation runs until the increase of the percentage likelihood falls below the `eps` value (default is `1e-9`), or when the maximum number of iterations is reached.
#' @param re_estimate Boolean. Indicates whether the model parameters should be re-estimated.
#' The default is `TRUE`. It can be set to `FALSE` to keep the model frozen for some period of time, although prolonged use of a frozen model is not recommended.
#'
#' @return An object of class `"JD3_DFMESTIMATES"` is returned. The following are
#'   returned invisibly as a list:
#' * `dfm` `[[1]]` an object of class `"JD3_DFMMODEL"` containing the estimated model parameters;
#' * `data` `[[2]]` the value of the `data` argument;
#' * `is_standardized` `[[3]]` the value of the `standardized` argument;
#' * `input_standardization` `[[4]]` the value of the `input_standardization` argument;
#' * `log_likelihood` `[[5]]` the estimated log-likelihood;
#' * `gradient` `[[6]]` the estimated gradient at the solution;
#' * `hessian` `[[7]]` the estimated hessian matrix at the solution;
#' * `has_converged` `[[8]]` a boolean indicating whether the ML estimation process has converged.
#'
#' @export
#'
#' @seealso `create_model()` to define a new model,
#'
#' `estimate_pca()` for estimation using principal components analysis,
#'
#' `estimate_em()` for estimation using EM algorithm,
#'
#' `get_results()` to access estimation results,
#'
#' `get_forecasts()` to obtain forecasts.
#'
#'
#' For more information, see the vignette:
#'
#' `utils::browseVignettes()`, e.g. `browseVignettes(package = "rjd3nowcasting")`
#'
#' @examplesIf rjd3jars::check_java_version(silent = TRUE)
#' # input data
#' set.seed(100)
#' data <- ts(matrix(rnorm(500), 100, 5),
#'            frequency = 12,
#'            start = c(2010, 1))
#' data[100, 1] <- data[99:100, 2] <- data[(1:100)[-seq(3, 100, 3)], 5] <- NA
#'
#' # define a new model
#' dfm <- create_model(
#'     nfactors = 2,
#'     nlags = 2,
#'     factors_type = c("M", "M", "YoY", "M", "Q"),
#'     factors_loading = matrix(data = TRUE, 5, 2),
#'     var_init = "Unconditional"
#' )
#'
#' # estimate by maximum likelihood
#' est_ml_a <- estimate_ml(dfm, data)
#'
#' # results and forecasts
#' rslts_a <- get_results(est_ml_a)
#' fcsts_a <- get_forecasts(est_ml_a)
#'
#' # no re-estimation of a previous model while integrating updated data in the output object
#' data_new <- data
#' data_new[99, 2] <- 1
#' est_ml_b <- estimate_ml(est_ml_a$dfm,
#'                         data_new,
#'                         input_standardization = rslts_a$preprocessing$sample_mean_stdev,
#'                         re_estimate = FALSE)
#'
estimate_ml <- function(dfm,
                        data,
                        standardized = FALSE,
                        input_standardization = NULL,
                        pca_init = TRUE,
                        em_init = TRUE,
                        em_max_iter = 100,
                        em_eps = 1e-9,
                        max_iter = 1000,
                        max_block_iter = 5,
                        simpl_model_iter = 15,
                        independent_var_shocks = FALSE,
                        mixedEstimation = TRUE,
                        eps = 1e-9,
                        re_estimate = TRUE) {

    if (re_estimate) {
        jdfm <- .r2jd_dfm(dfm)
        freq <- stats::frequency(data)
        start <- stats::start(data)
        jdata <- rjd3toolkit::.r2jd_matrix(data)
        if (is.null(input_standardization)) {
            standardization_mean <- standardization_stdev <- .jnull(class = "[D")
        } else {
            standardization_mean <- input_standardization[, 1]
            standardization_stdev <- input_standardization[, 2]
        }

        jest <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DfmEstimates;",
            "estimate_ML",
            jdfm$internal,
            jdata,
            as.integer(freq),
            .jarray(as.integer(start)),
            standardized,
            standardization_mean,
            standardization_stdev,
            pca_init,
            em_init,
            as.integer(em_max_iter),
            as.numeric(em_eps),
            as.integer(max_iter),
            as.integer(max_block_iter),
            as.integer(simpl_model_iter),
            independent_var_shocks,
            mixedEstimation,
            as.numeric(eps)
        )

        jmodel <- .jcall(
            "jdplus/dfm/base/r/DynamicFactorModels",
            "Ljdplus/dfm/base/core/DynamicFactorModel;",
            "getDfm",
            jest
        )

        dfm_list <- .jd2r_dfm(jmodel)
        jest <- rjd3toolkit::.jd3_object(jest, result = TRUE)
        ll <- rjd3toolkit::result(jest, "likelihood_ll")
        gradient <- rjd3toolkit::result(jest, "gradient")
        hessian <- rjd3toolkit::result(jest, "hessian")
        has_converged <- rjd3toolkit::result(jest, "has_converged")

    } else {
        if (!standardized && is.null(input_standardization)) {
            stop(
                "Re-estimation is turned off. Please provide the original mean and standard deviation in 'input_standardization'."
            )
        }

        dfm_list <- dfm
        ll <- NA
        gradient <- numeric()
        hessian <- matrix()
        has_converged <- "not re-estimated"
    }

    return(structure(
        list(
            dfm = dfm_list,
            data = data,
            is_standardized = standardized,
            input_standardization = input_standardization,
            log_likelihood = ll,
            gradient = gradient,
            hessian = hessian,
            has_converged = has_converged
        ),
        class = "JD3_DFMESTIMATES"
    ))
}


