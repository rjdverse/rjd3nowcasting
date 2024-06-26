% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{get_forecasts}
\alias{get_forecasts}
\title{Get DFM forecasts}
\usage{
get_forecasts(dfm_estimates, n_fcst = 3)
}
\arguments{
\item{dfm_estimates}{an object of class 'JD3_DfmEstimates'}

\item{n_fcst}{Integer. Number of forecast periods required.}
}
\value{
an object of class 'JD3_DfmForecasts'
}
\description{
Get DFM forecasts
}
\examples{
set.seed(100)
data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
dfm <- create_model(nfactors=2,
                    nlags=2,
                    factors_type = c("M", "M", "YoY", "M", "Q"),
                    factors_loading = matrix(data=TRUE, 5, 2),
                    var_init = "Unconditional")
est_em<-estimate_em(dfm, data)
fcst<-get_forecasts(est_em, n_fcst = 2)

}
