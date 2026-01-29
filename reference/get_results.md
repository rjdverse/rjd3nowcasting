# Get DFM results

Get DFM results

## Usage

``` r
get_results(dfm_estimates)
```

## Arguments

- dfm_estimates:

  an object of class 'JD3_DfmEstimates'

## Value

an object of class 'JD3_DfmResults'

## Examples

``` r
set.seed(100)
data<-ts(matrix(rnorm(500), 100, 5), frequency = 12, start = c(2010,1))
data[100,1]<-data[99:100,2]<-data[(1:100)[-seq(3,100,3)],5]<-NA
dfm <- create_model(nfactors=2,
                    nlags=2,
                    factors_type = c("M", "M", "YoY", "M", "Q"),
                    factors_loading = matrix(data=TRUE, 5, 2),
                    var_init = "Unconditional")
est_em<-estimate_em(dfm, data)
rslt_em<-get_results(est_em)
```
