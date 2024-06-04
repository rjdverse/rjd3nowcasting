#' @include estimation.R
DFMRESULTS<-'JD3_DfmResults'
DFMFORECASTS<-'JD3_DfmForecasts'

.dfmProcess<-function(dfm, data, standardized = FALSE, input_standardization = NULL){

  freq<-stats::frequency(data)
  start<-start(data)
  jdata<-rjd3toolkit::.r2jd_matrix(data)
  if(is.null(input_standardization)){
    standardization_mean <- standardization_stdev <- .jnull(class = "[D")
  } else{
    standardization_mean <- input_standardization[,1]
    standardization_stdev <- input_standardization[,2]
  }

  jrslts<-rjd3toolkit::.jd3_object(.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                                   "Ljdplus/dfm/base/core/DfmResults;",
                                   "process",
                                   .r2jd_dfm(dfm)$internal,
                                   jdata,
                                   as.integer(freq),
                                   .jarray(as.integer(start)),
                                   standardized,
                                   standardization_mean,
                                   standardization_stdev,
                                   as.integer(0)), result = TRUE)

  return(
    list(
      series_names=colnames(data),
      start=start,
      freq=freq,
      original_data=rjd3toolkit::result(jrslts, "input"),
      sample_mean=rjd3toolkit::result(jrslts, "sample_mean"),
      sample_stddev=rjd3toolkit::result(jrslts, "sample_stddev"),
      transformed_data=rjd3toolkit::result(jrslts, "input_transformed"),
      factors=rjd3toolkit::result(jrslts, "factors"),
      factors_stderr=rjd3toolkit::result(jrslts, "factors_stderr"),
      factors=rjd3toolkit::result(jrslts, "factors"),
      residuals=rjd3toolkit::result(jrslts, "residuals"),
      standard_residuals=rjd3toolkit::result(jrslts, "residuals_standardized")
    )
  )
}

# n_out = Number of out-of-periods considered in the processing. Must be
# greater than the required number of forecasts asked in the forecasting
# function.
#
.dfmProcessExt<-function(dfm, data, standardized = FALSE, input_standardization = NULL,
                         n_out = 12){

  freq <- stats::frequency(data)
  start <- start(data)
  jdata <- rjd3toolkit::.r2jd_matrix(data)
  if(is.null(input_standardization)){
    standardization_mean <- standardization_stdev <- .jnull(class = "[D")
  } else{
    standardization_mean <- input_standardization[,1]
    standardization_stdev <- input_standardization[,2]
  }
  jrslts<-rjd3toolkit::.jd3_object(.jcall("jdplus/dfm/base/r/DynamicFactorModels",
                                          "Ljdplus/dfm/base/core/DfmResults;",
                                          "process",
                                          .r2jd_dfm(dfm)$internal,
                                          jdata,
                                          as.integer(freq),
                                          .jarray(as.integer(start)),
                                          standardized,
                                          standardization_mean,
                                          standardization_stdev,
                                          as.integer(n_out)), result = TRUE)

  return(
    list(
      series_names=colnames(data),
      start=start,
      freq=freq,
      original_data=rjd3toolkit::result(jrslts, "input"),
      sample_mean=rjd3toolkit::result(jrslts, "sample_mean"),
      sample_stddev=rjd3toolkit::result(jrslts, "sample_stddev"),
      transformed_data=rjd3toolkit::result(jrslts, "input_transformed"),
      forecasts_T=rjd3toolkit::result(jrslts, paste0("forecasts_transformed(", n_out,")")),
      forecasts_T_stderr=rjd3toolkit::result(jrslts, paste0("forecasts_transformed_stderr(", n_out,")")),
      forecasts=rjd3toolkit::result(jrslts, paste0("forecasts(", n_out,")")),
      forecasts_stderr=rjd3toolkit::result(jrslts, paste0("forecasts_stderr(", n_out,")"))
    )
  )
}

.get_preprocessing <- function(dfm_rslts){

  series_names<-dfm_rslts$series_names

  # Original data
  data<-stats::ts(dfm_rslts$original_data, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  colnames(data)<-dfm_rslts$series_names

  # Sample mean and standard deviation
  sample_mean_stddev<-cbind(dfm_rslts$sample_mean, dfm_rslts$sample_stddev)
  colnames(sample_mean_stddev)<-c("sample_mean","sample_stddev")
  rownames(sample_mean_stddev)<-series_names

  # Transformed data
  data_t<-stats::ts(dfm_rslts$transformed_data, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  colnames(data_t)<-series_names

  return(list(
    original_data=data,
    sample_mean_stdev=sample_mean_stddev,
    transformed_data=data_t)
  )
}

.get_parameters <- function(dfm, series_names){

  # VAR coefficients
  var_coef<-dfm$var_coefficients
  nfactors<-nrow(var_coef)
  nlags<-ncol(var_coef)/nfactors

  var_coef_cnames<-var_coef_rnames<-character()
  k=1
  for(i in 1:nlags){
    for(j in 1:nfactors){
      if(i == 1) var_coef_rnames[j]<-paste0("F",j)
      var_coef_cnames[k]<-paste0("F",j,"[",-i,"]")
      k<-k+1
    }
  }
  colnames(var_coef)<-var_coef_cnames
  rownames(var_coef)<-var_coef_rnames

  # VAR var-cov of errors
  var_err_variance<-dfm$var_errors_variance
  colnames(var_err_variance)<-rownames(var_err_variance)<-var_coef_rnames

  # Measurement coefficients
  measurement_coef<-dfm$measurement_coefficients
  colnames(measurement_coef)<-var_coef_rnames
  rownames(measurement_coef)<-series_names

  # Variance of the idiosyncratic errors in measurement equation
  measurement_err_variance<-as.matrix(dfm$measurement_errors_variance, ncol = 1)
  colnames(measurement_err_variance)<-"idiosyncratic_variance"
  rownames(measurement_err_variance)<-series_names

  return(list(
    var_coefficients=var_coef,
    var_errors_variance=var_err_variance,
    measurement_coefficients=measurement_coef,
    measurement_errors_variance=measurement_err_variance)
  )
}

.get_factors <- function(dfm_rslts){

  # Factors
  factors<-stats::ts(dfm_rslts$factors, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  factors_cnames<-character()
  for(i in 1:ncol(factors)){
    factors_cnames[i]<-paste0("F",i)
  }
  colnames(factors)<-factors_cnames

  # Factors stdev
  factors_stdev<-stats::ts(dfm_rslts$factors_stderr, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  colnames(factors_stdev)<-factors_cnames

  return(list(
    factors=factors,
    factors_stdev=factors_stdev)
  )
}

# Residuals are the one-step forecast errors, while the standardized residuals
# are the one-step forecast errors after standardization (i.e. divided by
# their standard deviation). The basic diagnostics for normality,
# heteroscedasticity and auto-correlation should be performed on the
# standardized residuals.
.get_residuals <- function(dfm_rslts){
  res<-stats::ts(dfm_rslts$residuals, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  res_std<-stats::ts(dfm_rslts$standard_residuals, frequency=dfm_rslts$freq, start=dfm_rslts$start)
  colnames(res)<-colnames(res_std)<-dfm_rslts$series_names

  return(list(
    residuals=res,
    standardized_residuals=res_std)
  )
}

.get_likelihood <- function(dfm_estimates){

  ll<-dfm_estimates$log_likelihood
  has_converged<-dfm_estimates$has_converged

  return(list(
    log_likelihood=ll,
    has_converged=has_converged
  )
  )
}

#' Get DFM results
#' @param dfm_estimates an object of class 'JD3_DfmEstimates'
#' @return an object of class 'JD3_DfmResults'
#' @export
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
#' rslt_em<-get_results(est_em)
#'
get_results <- function(dfm_estimates){

  data<-dfm_estimates$data
  dfm<-dfm_estimates$dfm
  is_standardized<-dfm_estimates$is_standardized
  input_standardization<-dfm_estimates$input_standardization
  dfm_rslts<-.dfmProcess(dfm, data, is_standardized, input_standardization)

  return(structure(list(
    preprocessing=.get_preprocessing(dfm_rslts),
    parameters=.get_parameters(dfm, dfm_rslts$series_names),
    factors=.get_factors(dfm_rslts),
    residuals=.get_residuals(dfm_rslts),
    likelihood=.get_likelihood(dfm_estimates)),
    class = DFMRESULTS))
}

#' Print function for objects of class 'JD3_DfmResults'
#'
#' @param x an object of class 'JD3_DfmResults'
#' @param \dots further arguments passed to the print() function.
#' @export
#' @exportS3Method print JD3_DfmResults
#'
print.JD3_DfmResults <- function(x, ...){

  sample_mean_stdev<-round(x$preprocessing$sample_mean_stdev,3)
  param_mcoeff<-round(x$parameters$measurement_coefficients,3)
  param_mvar<-round(x$parameters$measurement_errors_variance,3)
  loadings<-cbind(sample_mean_stdev, param_mcoeff, param_mvar)
  colnames(loadings)<-c("Sample mean", "Stdev",
                        paste0("Coeff. ", colnames(param_mcoeff)), "Idiosyncratic variance")

  print(list(loadings=loadings,
             VAR_model=round(x$parameters$var_coefficients,3),
             Innovative_variance=round(x$parameters$var_errors_variance,3), ...))
}



#' Get DFM forecasts
#'
#' @param dfm_estimates an object of class 'JD3_DfmEstimates'
#' @param n_fcst Integer. Number of forecast periods required.
#' @return an object of class 'JD3_DfmForecasts'
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
#' fcst<-get_forecasts(est_em, n_fcst = 2)
#'
get_forecasts <- function(dfm_estimates, n_fcst = 3){

  if(n_fcst < 1) n_fcst<-1
  data<-dfm_estimates$data
  dfm<-dfm_estimates$dfm
  is_standardized<-dfm_estimates$is_standardized
  input_standardization<-dfm_estimates$input_standardization
  dfm_rslts<-.dfmProcessExt(dfm, data, is_standardized,
                            input_standardization, n_fcst)

  # Transformed series
  freq<-dfm_rslts$freq
  start<-dfm_rslts$start
  fcsts_t<-stats::ts(dfm_rslts$forecasts_T, frequency=freq, start=start)
  fcsts_t_stderr<-stats::ts(dfm_rslts$forecasts_T_stderr, frequency=freq, start=start)
  colnames(fcsts_t)<-colnames(fcsts_t_stderr)<-dfm_rslts$series_names

  # Original series
  fcsts<-stats::ts(dfm_rslts$forecasts, frequency=freq, start=start)
  fcsts_stderr<-stats::ts(dfm_rslts$forecasts_stderr, frequency=freq, start=start)
  colnames(fcsts)<-colnames(fcsts_stderr)<-dfm_rslts$series_names

  # Restrict output to forecasts only
  data<-stats::ts(dfm_rslts$original_data, frequency=freq, start=start)

  nc<-ncol(data)
  nf<-vector(mode="integer", length=nc)
  for(j in 1:nc){
    cj<-data[,j]
    nf[j]<-min(which(!is.na(rev(cj))))-1
  }

  nf_max<-max(nf)
  strt<-stats::time(fcsts)[nrow(data)-nf_max+1]
  strt_yr<-floor(strt)
  strt_mth<-round((strt%%1)*freq+1,0)

  fcsts_only<-stats::window(fcsts, start = c(strt_yr, strt_mth))
  fcsts_only_stderr<-stats::window(fcsts_stderr, start = c(strt_yr, strt_mth))

  if(nf_max > 0){
    for(j in 1:nc){
      n_na<-nf_max-nf[j]
      if(n_na>0) fcsts_only[1:n_na,j]<-fcsts_only_stderr[1:n_na,j]<-NA
    }
  }

  return(structure(list(
    transformed_forecasts=fcsts_t,
    transformed_forecasts_stdev=fcsts_t_stderr,
    forecasts=fcsts,
    forecasts_stdev=fcsts_stderr,
    forecasts_only=fcsts_only,
    forecasts_only_stdev=fcsts_only_stderr),
    class = DFMFORECASTS))
}

#' Print function for objects of class 'JD3_DfmForecasts'
#'
#' @param x an object of class 'JD3_DfmForecasts'
#' @param \dots further arguments passed to the print() function.
#' @export
#' @exportS3Method print JD3_DfmForecasts
#'
print.JD3_DfmForecasts <- function(x, ...){
  print(list(forecasts_only=x$forecasts_only, ...))
}

#' Plot function for objects of class 'JD3_DfmForecasts'
#'
#' @param x an object of class 'JD3_DfmForecasts'
#' @param series_name Character. Name of the series to plot. By default, the
#'   first series will be plotted.
#' @param \dots further arguments passed to ts.plot().
#' @exportS3Method plot JD3_DfmForecasts
#' @export
#'
plot.JD3_DfmForecasts <- function(x, series_name = NULL, ...){

  fcst<-x$forecasts
  fcst_stdev<-x$forecasts_stdev
  fcst_only<-x$forecasts_only

  if(is.null(series_name)){
    series_name<-colnames(fcst)[1]
  }

  if(series_name %in% colnames(fcst)){
    s<-fcst[,series_name]
    s_lb<-s - 1.28 * fcst_stdev[,series_name]
    s_ub<-s + 1.28 * fcst_stdev[,series_name]
    sf<-fcst_only[,series_name]
  }else{
    stop("series name not found!")
  }

  stats::ts.plot(s_lb, s_ub, s, sf,
                 gpars=list(main=series_name, sub="Forecasts with a 80% PI", xlab="", ylab= "", lty=c(3, 3, 1, 1), xaxt="n", type = "o", pch=20, cex=0.8, las=2, col = c("orange","orange","black","red"), ...))
  graphics::axis(1, at=seq(stats::start(s)[1], stats::end(s)[1], by = 1), las=2)
  graphics::legend("topleft", legend=c("series", "forecasts", "80% PI"), col=c("black","red","orange"), lty=c(1, 1, 3), cex=0.8)
}




