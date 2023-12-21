## Overview

Nowcasting is often defined as the prediction of the present, the very near future and the very recent past. 

rjd3nowcasting provides helps to operationalize the process of nowcasting. This first version can be used to specify and estimate dynamic factor models. 
A later version is expected to include the concept of "news" similar to the [Nowcasting plugin](https://github.com/nbbrd/jdemetra-nowcasting/tree/master) 
of the Graphical User Interface of JDemetra+ v2.

## Installation

To get the current stable version (from the latest release):

``` r
# install.packages("remotes")
remotes::install_github("rjdemetra/rjd3toolkit@*release")
remotes::install_github("rjdemetra/rjd3nowcasting@*release")
```

To get the current development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("rjdemetra/rjd3nowcasting")
```

## Usage

### Input
``` r
set.seed(100)
data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
```

### Model
``` r
dfm_model <- rjd3nowcasting::model(nfactors=2,
                                   nlags=2,
                                   factors_type = c("M", "M", "YoY", "M", "Q"),
                                   factors_loading = matrix(data=TRUE, 5, 2),
                                   var_init = "Unconditional")
```

### Estimation
``` r
rslt_ml<-rjd3nowcasting::estimate_ml(dfm_model, data)
# or rslt_em<-rjd3nowcasting::estimate_em(dfm_model, data)
# or rslt_pca<-rjd3nowcasting::estimate_pca(dfm_model, data)
```

### Results
``` r
fcst<-rjd3nowcasting::get_forecasts(rslt_ml,nf = 2,forecasts_only = TRUE)
params<-rjd3nowcasting::get_parameters(rslt_ml)
factors<-rjd3nowcasting::get_factors(rslt_ml)
#...

print(rslt_ml)
summary(rslt_ml)
plot(rslt_ml)
```

