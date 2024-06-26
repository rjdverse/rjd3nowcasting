% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/news.R
\name{get_news}
\alias{get_news}
\title{DFM News analysis}
\usage{
get_news(dfm_estimates, new_data, target_series = NULL, n_fcst = 3)
}
\arguments{
\item{dfm_estimates}{an object of class 'JD3_DfmEstimates'. Typically
generated by the functions estimate_pca(), estimate_em() or estimate_ml().}

\item{new_data}{an mts object containing the updated dataset.}

\item{target_series}{the name of the series of interest. By default, the
first series is considered.}

\item{n_fcst}{the number of forecasting periods to consider. Default is 3.}
}
\value{
An object of class 'JD3_DfmNews'
}
\description{
DFM News analysis
}
\examples{
set.seed(100)
data_t1<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
data_t1[100,1]<-data_t1[99:100,2]<-data_t1[(1:100)[-seq(3,100,3)],5]<-NA
data_t2<-ts(rbind(data_t1, rep(NA,5)), frequency = 12, start = c(2010,1))
data_t2[100,1]<-data_t2[99,2]<-data_t2[101,3]<-data_t2[101,4]<-1

dfm_model <- create_model(nfactors=2,
                          nlags=2,
                          factors_type = c("M", "M", "YoY", "M", "Q"),
                          factors_loading = matrix(TRUE, 5, 2),
                          var_init = "Unconditional")

est_em<-estimate_em(dfm_model, data_t1)
# or to use any previous frozen model:
# est_em_frozen<-estimate_em(dfm_model, data_t1, re_estimate = FALSE)

news<-get_news(est_em, data_t2, target_series = "Series 2", n_fcst = 2)

}
\references{
Banbura and Modugno (2010) - Maximum likelihood estimation of factor models on data sets with arbitrary pattern of missing data
}
