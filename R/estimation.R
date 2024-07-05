DFMMODEL<-'JD3_DfmModel'
DFMESTIMATES<-'JD3_DfmEstimates'


.r2jd_dfm<-function(dfm_list){
  # pm dfm_list = list(var_coef, var_var, m_coef, m_var, factors_type, var_init)

  m_coef <- dfm_list$measurement_coefficients
  m_var <- dfm_list$measurement_errors_variance
  var_coef <- dfm_list$var_coefficients
  var_var <- dfm_list$var_errors_variance
  factors_type <- dfm_list$factors_type
  var_init <- dfm_list$initialization_type

  nf<-ncol(m_coef)
  nl<-ncol(var_coef)/nf
  loadings<-ifelse(is.nan(m_coef), FALSE, TRUE)
  jloadings <- rjd3toolkit::.r2jd_matrix(loadings)
  jvar_coef <- rjd3toolkit::.r2jd_matrix(var_coef)
  jvar_var <- rjd3toolkit::.r2jd_matrix(var_var)
  jm_coef <- rjd3toolkit::.r2jd_matrix(m_coef)

  jmodel<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
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
                 .jarray(as.numeric(m_var)))

  return(rjd3toolkit::.jd3_object(jmodel, DFMMODEL, result = TRUE))
}

.jd2r_dfm<-function(jmodel){

  jmodel<-rjd3toolkit::.jd3_object(jmodel, result = TRUE)

  vc<-rjd3toolkit::result(jmodel,"var_coefficients")
  vv<-rjd3toolkit::result(jmodel,"var_errors_variance")
  mc<-rjd3toolkit::result(jmodel,"measurement_coefficients")
  mv<-rjd3toolkit::result(jmodel,"measurement_errors_variance")
  init<-rjd3toolkit::result(jmodel,"initialization_type")
  ftypes<-rjd3toolkit::result(jmodel,"factors_type") # number of lags implied by the measurement
  ftypes_char<-sapply(as.character(ftypes), switch, "1"="M", "5"="Q", "12"="YoY",
                      USE.NAMES = FALSE)

  return(
    structure(list(
      var_coefficients=vc,
      var_errors_variance=vv,
      measurement_coefficients=mc,
      measurement_errors_variance=mv,
      initialization_type=init,
      factors_type=ftypes_char),
      class = DFMMODEL)
  )
}


#' Create Dynamic Factor Model
#'
#' @param nfactors Integer. Number of factors.
#' @param nlags Integer. Number of lags in VAR equations.
#' @param factors_type Character vector. Respecting the order of the series in
#'   the input data, you must refer here the link between the (transformed)
#'   series and the factors. Three options are possible:
#'  \itemize{
#'    \item{"M"}{: Variables expressed in terms of monthly growth rates can be
#'    linked to a factor representing the underlying monthly growth rate of the
#'    economy if "M" is selected}
#'    \item{"Q"}{: Monthly or quarterly variables that are correlated with the the
#'    underlying quarterly growth rate of the economy can be linked to a
#'    weighted average of the factors representing the underlying monthly
#'    growth rate of the economy. Such a weighted average is meant to represent
#'    quarterly growth rates, and it is implemented by selecting "Q"}
#'    \item{"YoY"}{: The variables can also be linked to the cumulative sum of the
#'    last 12 monthly factors. If the model is designed in such a way that the
#'    monthly factors represent monthly growth rates, the resulting cumulative
#'    sum boils down to the year-on-year growth rate. Thus, variables expressed
#'    in terms of year-on-year growth rates or surveys that are correlated with
#'    the year-on-year growth rates of the reference series should be linked to
#'    the factors using "YoY".}
#'  }
#' @param factors_loading Boolean matrix. It represents the factor loading
#'   structure. The dimension of the matrix should be 'number of series' x
#'   'number of factors'. For each row representing each series, the user must
#'   mention whether the corresponding factor loads on this series.
#' @param var_init Character. The first unobserved factors values in the sample
#'   is assumed to be either equal to zero or consistent with a normal
#'   distribution with mean zero and a variance corresponding to the
#'   unconditional variance of the VAR. The latter is the default.
#' @param var_coefficients Matrix. The default is NULL meaning that the VAR
#'   coefficients will be estimated from scratch. Alternatively, a matrix of
#'   pre-defined values can be passed in. Those would come typically from a
#'   previous model estimate and will serve as a starting point for the
#'   estimation step. The format of the matrix should be the same as the one
#'   produced by default by the create_model() function while keeping the
#' `var_coefficients` argument to its default value NULL.
#' @param var_errors_variance Matrix. The default is NULL meaning that the VAR
#'   errors variance will be estimated from scratch. Alternatively, a matrix of
#'   pre-defined values can be passed in. Those would come typically from a
#'   previous model estimate and will serve as a starting point for the
#'   estimation step. The format of the matrix should be the same as the one
#'   produced by default by the create_model() function while keeping the
#' `var_errors_variance` argument to its default value NULL.
#' @param measurement_coefficients Matrix. The default is NULL meaning that the
#'   measurement coefficients will be estimated from scratch. Alternatively, a
#'   matrix of pre-defined values can be passed in. Those would come typically
#'   from a previous model estimate and will serve as a starting point for the
#'   estimation step. The format of the matrix should be the same as the one
#'   produced by default by the create_model() function while keeping the
#' `measurement_coefficients` argument to its default value NULL.
#' @param measurement_errors_variance Numeric vector. The default is NULL
#'   meaning that the measurement errors variance will be estimated from
#'   scratch. Alternatively, a vector of pre-defined values can be passed in.
#'   Those would come typically from a previous model estimate and will serve as
#'   a starting point for the estimation step. The format of the vector should
#'   be the same as the one produced by default by the create_model() function
#'   while keeping the
#' `measurement_errors_variance` argument to its default value NULL.
#' @return an object of class 'JD3_DfmModel'
#' @export
#'
#' @examples
#'
#' # From scratch
#' dfm1 <- create_model(nfactors=2,
#'                      nlags=2,
#'                      factors_type = c("M", "M", "YoY", "M", "Q"),
#'                      factors_loading = matrix(data=TRUE, 5, 2),
#'                      var_init = "Unconditional")
#'
#' # From a previous estimate
#' set.seed(100)
#' data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
#' data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
#' est1<-estimate_em(dfm1, data)
#'
#' dfm2 <- create_model(nfactors=2,
#'                      nlags=2,
#'                      factors_type = c("M", "M", "YoY", "M", "Q"),
#'                      factors_loading = matrix(data=TRUE, 5, 2),
#'                      var_init = "Unconditional",
#'                      var_coefficients = est1$dfm$var_coefficients,
#'                      var_errors_variance = est1$dfm$var_errors_variance,
#'                      measurement_coefficients = est1$dfm$measurement_coefficients,
#'                      measurement_errors_variance = est1$dfm$measurement_errors_variance)
#' #est2<-estimate_em(dfm2, data)
#'
create_model<-function(nfactors, nlags, factors_type, factors_loading, var_init = c("Unconditional", "Zero"),
                       var_coefficients = NULL, var_errors_variance = NULL,
                       measurement_coefficients = NULL, measurement_errors_variance = NULL){

  var_init<-match.arg(var_init)
  jfactors_loading <- rjd3toolkit::.r2jd_matrix(factors_loading)
  jvar_coefficients <- rjd3toolkit::.r2jd_matrix(var_coefficients)
  jvar_errors_variance <- rjd3toolkit::.r2jd_matrix(var_errors_variance)
  jmeasurement_coefficients <- rjd3toolkit::.r2jd_matrix(measurement_coefficients)
  if (is.null(measurement_errors_variance)){
    jmeasurement_errors_variance <- .jnull("[D")
  } else {
    jmeasurement_errors_variance <- .jarray(as.numeric(measurement_errors_variance))
  }

  jmodel<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                 "Ljdplus/dfm/base/core/DynamicFactorModel;",
                 "model",
                 as.integer(nfactors),
                 as.integer(nlags),
                 .jarray(as.character(factors_type)),
                 jfactors_loading,
                 .jnew("java/lang/String", as.character(var_init)),
                 jvar_coefficients,
                 jvar_errors_variance,
                 jmeasurement_coefficients,
                 jmeasurement_errors_variance)

  return(.jd2r_dfm(jmodel))
}


#' Estimate DFM with Principal components Analysis
#'
#' @param dfm an object of class 'JD3_DfmModel'. Typically generated by the
#'   create_model() function.
#' @param data an mts object.
#' @param standardized Boolean. Indicate whether the input series were already
#'   standardized or not. Default is FALSE, meaning that a standardization of
#'   the series will be preliminary applied as part of the process.
#' @param input_standardization Matrix. Mean and standard deviation of the
#'   variables to consider for the pre-processing step of standardization.
#'   Default is NULL, meaning that they will be re-calculated based on the data.
#'   Typically, it can be filled with the output of the function
#'   `get_results()$preprocessing$sample_mean_stdev` applied on a
#'   previous estimate of the model. If provided manually, it must be a two
#'   columns matrix with the mean in the first column and the standard deviation
#'   in the second column. In the rows, the order of the variables should also
#'   be respected (similar to the data). Note that this argument must be filled
#'   if the re_estimate argument is set to FALSE. On the other hand, it is
#'   ignored if the standardized argument is set to TRUE.
#' @param re_estimate Boolean. Indicate whether the model will be re-estimated
#'   or not. Default is TRUE. Could be set to FALSE if, for some reasons during
#'   the production process, we wanted to freeze to model for some periods of
#'   time. It is not recommended to freeze the model for a long period.
#' @return an object of class 'JD3_DfmEstimates'
#' @export
#'
#' @examples
#' set.seed(100)
#' data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
#' data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
#' dfm <- create_model(nfactors=2,
#'                     nlags=2,
#'                     factors_type = c("M", "M", "YoY", "M", "Q"),
#'                     factors_loading = matrix(data=TRUE, 5, 2),
#'                     var_init = "Unconditional")
#' est_pca<-estimate_pca(dfm, data)
#'
#' #est_pca<-estimate_pca(dfm, data, re_estimate=FALSE) # model not re-estimated
#'
estimate_pca<-function(dfm, data, standardized = FALSE,
                       input_standardization = NULL, re_estimate = TRUE){

  if (re_estimate){
    jdfm<-.r2jd_dfm(dfm)
    freq<-stats::frequency(data)
    start<-start(data)
    jdata<-rjd3toolkit::.r2jd_matrix(data)
    if (is.null(input_standardization)){
      standardization_mean <- standardization_stdev <- .jnull(class = "[D")
    } else {
      standardization_mean <- input_standardization[,1]
      standardization_stdev <- input_standardization[,2]
    }

    jest<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                 "Ljdplus/dfm/base/core/DfmEstimates;",
                 "estimate_PCA",
                 jdfm$internal,
                 jdata,
                 as.integer(freq),
                 .jarray(as.integer(start)),
                 standardized,
                 standardization_mean,
                 standardization_stdev)

    jmodel<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                   "Ljdplus/dfm/base/core/DynamicFactorModel;",
                   "getDfm",
                   jest)

    dfm_list<-.jd2r_dfm(jmodel)
    jest<-rjd3toolkit::.jd3_object(jest, result = TRUE)
    ll<-rjd3toolkit::result(jest,"likelihood_ll")
    gradient<-rjd3toolkit::result(jest,"gradient")
    hessian<-rjd3toolkit::result(jest,"hessian")
    has_converged<-rjd3toolkit::result(jest,"has_converged")

  } else {
    if (!standardized && is.null(input_standardization)){
      stop("Since you chose not to re-estimate your model, you must also
           fill the argument 'input_standardization' with the original mean
           and standard deviation that was previously used to standardize your
           data beforehand")
    }

    dfm_list<-dfm
    ll<-NA
    gradient<-numeric()
    hessian<-matrix()
    has_converged<-"not re-estimated"
  }

  return(
    structure(list(
      dfm=dfm_list,
      data=data,
      is_standardized=standardized,
      input_standardization=input_standardization,
      log_likelihood=ll,
      gradient=gradient,
      hessian=hessian,
      has_converged=has_converged),
      class = DFMESTIMATES)
  )
}

#' Estimate DFM with Expectations-Maximization algorithm
#'
#' @param dfm an object of class 'JD3_DfmModel'. Typically generated by the
#'   create_model() function.
#' @param data an mts object.
#' @param standardized Boolean. Indicate whether the input series were already
#'   standardized or not. Default is FALSE, meaning that a standardization of
#'   the series will be preliminary applied as part of the process.
#' @param input_standardization Matrix. Mean and standard deviation of the
#'   variables to consider for the pre-processing step of standardization.
#'   Default is NULL, meaning that they will be re-calculated based on the data.
#'   Typically, it can be filled with the output of the function
#'   `get_results()$preprocessing$sample_mean_stdev` applied on a
#'   previous estimate of the model. If provided manually, it must be a two
#'   columns matrix with the mean in the first column and the standard deviation
#'   in the second column. In the rows, the order of the variables should also
#'   be respected (similar to the data). Note that this argument must be filled
#'   if the re_estimate argument is set to FALSE. On the other hand, it is
#'   ignored if the standardized argument is set to TRUE.
#' @param pca_init Boolean. Indicate whether a principal components analysis
#'   is performed beforehand and used as initial condition for the EM
#'   algorithm.
#' @param max_iter Integer. Maximum number of iterations.
#' @param eps Numeric. EM algorithm is run until the percentage likelihood does
#'   not increase by more than the eps value (1e-9 is the default) or until the
#'   maximum number of iterations is hit.
#' @param re_estimate Boolean. Indicate whether the model will be re-estimated
#'   or not. Default is TRUE. Could be set to FALSE if, for some reasons during
#'   the production process, we wanted to freeze to model for some periods of
#'   time. It is not recommended to freeze the model for a long period.
#' @return an object of class 'JD3_DfmEstimates'
#' @export
#'
#' @examples
#' set.seed(100)
#' data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
#' data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
#' dfm <- create_model(nfactors=2,
#'                     nlags=2,
#'                     factors_type = c("M", "M", "YoY", "M", "Q"),
#'                     factors_loading = matrix(data=TRUE, 5, 2),
#'                     var_init = "Unconditional")
#' est_em<-estimate_em(dfm, data)
#'
#' #est_em<-estimate_em(dfm, data, re_estimate=FALSE) # model not re-estimated
#'
estimate_em<-function(dfm, data, standardized = FALSE, input_standardization = NULL,
                      pca_init = TRUE, max_iter = 100, eps = 1e-9,
                      re_estimate = TRUE){

  if (re_estimate){
    jdfm<-.r2jd_dfm(dfm)
    freq<-stats::frequency(data)
    start<-start(data)
    jdata<-rjd3toolkit::.r2jd_matrix(data)
    if (is.null(input_standardization)){
      standardization_mean <- standardization_stdev <- .jnull(class = "[D")
    } else {
      standardization_mean <- input_standardization[,1]
      standardization_stdev <- input_standardization[,2]
    }

    jest<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
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
                 as.numeric(eps))

    jmodel<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                   "Ljdplus/dfm/base/core/DynamicFactorModel;",
                   "getDfm",
                   jest)

    dfm_list<-.jd2r_dfm(jmodel)
    jest<-rjd3toolkit::.jd3_object(jest, result = TRUE)
    ll<-rjd3toolkit::result(jest,"likelihood_ll")
    gradient<-rjd3toolkit::result(jest,"gradient")
    hessian<-rjd3toolkit::result(jest,"hessian")
    has_converged<-rjd3toolkit::result(jest,"has_converged")

  } else {
    if (!standardized && is.null(input_standardization)){
      stop("Since you chose not to re-estimate your model, you must also
           fill the argument 'input_standardization' with the original mean
           and standard deviation that was previously used to standardize your
           data beforehand")
    }

    dfm_list<-dfm
    ll<-NA
    gradient<-numeric()
    hessian<-matrix()
    has_converged<-"not re-estimated"
  }

  return(
    structure(list(
      dfm=dfm_list,
      data=data,
      is_standardized=standardized,
      input_standardization=input_standardization,
      log_likelihood=ll,
      gradient=gradient,
      hessian=hessian,
      has_converged=has_converged),
      class = DFMESTIMATES)
  )
}


#' Estimate DFM with Maximum Likelihood
#'
#' @param dfm an object of class 'JD3_DfmModel'. Typically generated by the
#'   create_model() function.
#' @param data an mts object.
#' @param standardized Boolean. Indicate whether the input series were already
#'   standardized or not. Default is FALSE, meaning that a standardization of
#'   the series will be preliminary applied as part of the process.
#' @param input_standardization Matrix. Mean and standard deviation of the
#'   variables to consider for the pre-processing step of standardization.
#'   Default is NULL, meaning that they will be re-calculated based on the data.
#'   Typically, it can be filled with the output of the function
#'   `get_results()$preprocessing$sample_mean_stdev` applied on a
#'   previous estimate of the model. If provided manually, it must be a two
#'   columns matrix with the mean in the first column and the standard deviation
#'   in the second column. In the rows, the order of the variables should also
#'   be respected (similar to the data). Note that this argument must be filled
#'   if the re_estimate argument is set to FALSE. On the other hand, it is
#'   ignored if the standardized argument is set to TRUE.
#' @param pca_init Boolean. Indicate whether a principal components analysis is
#'   performed beforehand and used as initial condition for either the EM
#'   algorithm (if em_init=TRUE) or directly for the ML estimation.
#' @param em_init Boolean. Indicate whether the EM algorithm is performed
#'   beforehand and used as initial condition for the ML estimation.
#' @param em_max_iter Integer. Maximum number of iterations of the EM algorithm.
#'   Ignored if em_init = FALSE.
#' @param em_eps Numeric. EM algorithm is run until the percentage likelihood
#'   does not increase by more than the eps value (1e-9 is the default) or until
#'   the maximum number of iterations is hit. Ignored if em_init = FALSE.
#' @param max_iter Integer. Maximum number of iterations for the ML estimation.
#' @param max_block_iter Integer. Maximum number of iterations in optimization
#'   by block. The model parameters are divided in two blocks: one related to
#'   the measurement equations and one to the VAR equations. While the EM
#'   algorithm requires one iteration per block, the numerical optimization
#'   allows us to set the number of iterations desired per block.
#' @param simpl_model_iter Integer. Number of simplified model iterations
#'   allowed.
#' @param independent_var_shocks Boolean. Whether we assume that shocks in the
#'   VAR block are independent.
#' @param mixedEstimation Boolean. The mixed estimation option alternates
#'   between the iterations for the VAR block alone and simultaneous iterations
#'   for the two blocks.
#' @param eps Numeric. ML estimation is run until the percentage likelihood does
#'   not increase by more than the eps value (1e-9 is the default) or until the
#'   maximum number of iterations is hit.
#' @param re_estimate Boolean. Indicate whether the model will be re-estimated
#'   or not. Default is TRUE. Could be set to FALSE if, for some reasons during
#'   the production process, we wanted to freeze to model for some periods of
#'   time. It is not recommended to freeze the model for a long period.
#' @return an object of class 'JD3_DfmEstimates'
#' @export
#'
#' @examples
#' set.seed(100)
#' data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
#' data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
#' dfm <- create_model(nfactors=2,
#'                     nlags=2,
#'                     factors_type = c("M", "M", "YoY", "M", "Q"),
#'                     factors_loading = matrix(data=TRUE, 5, 2),
#'                     var_init = "Unconditional")
#' est_ml<-estimate_ml(dfm, data)
#'
#' #est_ml<-estimate_ml(dfm, data, re_estimate=FALSE) # model not re-estimated
#'
estimate_ml<-function(dfm, data, standardized = FALSE, input_standardization = NULL,
                      pca_init = TRUE, em_init = TRUE, em_max_iter = 100,
                      em_eps = 1e-9, max_iter = 1000, max_block_iter = 5,
                      simpl_model_iter = 15, independent_var_shocks = FALSE,
                      mixedEstimation = TRUE, eps=1e-9, re_estimate = TRUE){

  if (re_estimate){
    jdfm<-.r2jd_dfm(dfm)
    freq<-stats::frequency(data)
    start<-start(data)
    jdata<-rjd3toolkit::.r2jd_matrix(data)
    if (is.null(input_standardization)){
      standardization_mean <- standardization_stdev <- .jnull(class = "[D")
    } else {
      standardization_mean <- input_standardization[,1]
      standardization_stdev <- input_standardization[,2]
    }

    jest<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
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
                 as.numeric(eps))

    jmodel<-.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                   "Ljdplus/dfm/base/core/DynamicFactorModel;",
                   "getDfm",
                   jest)

    dfm_list<-.jd2r_dfm(jmodel)
    jest<-rjd3toolkit::.jd3_object(jest, result = TRUE)
    ll<-rjd3toolkit::result(jest,"likelihood_ll")
    gradient<-rjd3toolkit::result(jest,"gradient")
    hessian<-rjd3toolkit::result(jest,"hessian")
    has_converged<-rjd3toolkit::result(jest,"has_converged")

  } else {
    if (!standardized && is.null(input_standardization)){
      stop("Since you chose not to re-estimate your model, you must also
           fill the argument 'input_standardization' with the original mean
           and standard deviation that was previously used to standardize your
           data beforehand")
    }

    dfm_list<-dfm
    ll<-NA
    gradient<-numeric()
    hessian<-matrix()
    has_converged<-"not re-estimated"
  }

  return(
    structure(list(
      dfm=dfm_list,
      data=data,
      is_standardized=standardized,
      input_standardization=input_standardization,
      log_likelihood=ll,
      gradient=gradient,
      hessian=hessian,
      has_converged=has_converged),
      class = DFMESTIMATES)
  )
}

#' Print function for objects of class 'JD3_DfmEstimates'
#'
#' @param x an object of class 'JD3_DfmEstimates'
#' @param \dots further arguments passed to the print() function.
#' @export
#' @exportS3Method print JD3_DfmEstimates
#'
print.JD3_DfmEstimates <- function(x, ...){

  print(list(has_converged=toupper(x$has_converged),
             log_likelihood=x$log_likelihood, ...))

}
