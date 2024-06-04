
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rjd3nowcasting

<!-- badges: start -->
<!-- badges: end -->

## Overview

Nowcasting is often defined as the prediction of the present, the very
near future and the very recent past.

The package *rjd3nowcasting* provides helps to operationalize the process of
nowcasting. This package can be used to specify and estimate Dynamic Factor Models. Latest version of the package also includes news analysis. This R package uses the efficient libraries of [JDemetra+ v3](https://github.com/jdemetra/jdplus-nowcasting). 

The way the package was conceived is inspired by the [GUI add-in](https://github.com/nbbrd/jdemetra-nowcasting) developed for JDemetra+ V2 and it provides about the same functionalities (except for the real-time simulation), but in the flexible R environment.

## Installation

Running rjd3 packages requires **Java 17 or higher**. How to set up such
a configuration in R is explained
[here](https://jdemetra-new-documentation.netlify.app/#Rconfig)

To get the current stable version (from the latest release):

``` r
# install.packages("remotes")
remotes::install_github("rjdverse/rjd3toolkit@*release")
remotes::install_github("rjdverse/rjd3nowcasting@*release")
```

To get the current development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("rjdverse/rjd3nowcasting")
```

## Usage

``` r
library("rjd3nowcasting")
```

Once the package is loaded, there are four steps to follow:

1. Import data 
2. Create or update the model 
3. Estimate the model 
4. Get results

### 1. Data

``` r
set.seed(100)
data0 <- stats::ts(
    data = matrix(rnorm(500), 100, 5), 
    frequency = 12, 
    start = c(2010, 1)
)
data0[100, 1] <- data0[99:100, 2] <- data0[(1:100)[-seq(3, 100, 3)], 5] <- NA

data1 <- stats::ts(
    data = rbind(data0, c(NA, NA, 1, 1, NA)), 
    frequency = 12, 
    start = c(2010, 1)
)
data1[100,1] <- data1[99,2] <- 1
```

### 2. Create or update the model

``` r
# Create model from scratch
dfm0 <- create_model(nfactors=2,
                     nlags=2,
                     factors_type = c("M", "M", "YoY", "M", "Q"),
                     factors_loading = matrix(data = TRUE, 5, 2),
                     var_init = "Unconditional")

# Update model
est0 <- estimate_em(dfm0, data0) # cfr. next step

dfm1 <- est0$dfm # R object (list) to potentially save from one time to another  
# or, equivalently,
dfm1 <- create_model(nfactors=2,
                     nlags=2,
                     factors_type = c("M", "M", "YoY", "M", "Q"),
                     factors_loading = matrix(data = TRUE, 5, 2),
                     var_init = "Unconditional",
                     var_coefficients = est0$dfm$var_coefficients,
                     var_errors_variance = est0$dfm$var_errors_variance,
                     measurement_coefficients = est0$dfm$measurement_coefficients,
                     measurement_errors_variance = est0$dfm$measurement_errors_variance)
```

### 3. Estimate the model

``` r
est1 <- estimate_ml(dfm1, data1)
# or est1<-estimate_em(dfm1, data1)
# or est1<-estimate_pca(dfm1, data1)
```

### 4. Get results

``` r
rslt1 <- get_results(est1)
print(rslt1)
fcst1 <- get_forecasts(est1, n_fcst = 2)
print(fcst1)
plot(fcst1)

news1 <- get_news(est0, data1, target_series = "Series 1", n_fcst = 2)
print(news1)
plot(news1)
```

## Package Maintenance and contributing

Any contribution is welcome and should be done through pull requests
and/or issues. pull requests should include **updated tests** and
**updated documentation**. If functionality is changed, docstrings
should be added or updated.

## Licensing

The code of this project is licensed under the [European Union Public
Licence
(EUPL)](https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12).
