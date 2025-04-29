# BIOS735 Group 2 R Package


## Installation

``` r
# Note: this will install the package 'pastasales'
# Note: Compilation of Rcpp may be required
remotes::install_github("jtlandis/bios735-group-proj")
```

## Final Report
The BIOS735 final report can be found in `report/` as Report.pdf with 
related Overleaf documentation in `report/overleaf/` folder. 


## pastasales

``` r
library(pastasales) # main package
library(dplyr) # helpful for analysis
```

``` r
# Example: Fitting a Poisson Autoregressive (PAR) Model
# The `pastasales` package provides tools for analyzing hierarchical sales data
# using Poisson autoregressive models. Below is a quick example of how to use
# the package to fit a PAR model.

# Step 1: Subset and prepare the data
model_spec <- data_set_tidy |>
  # Select a smaller subset for faster demonstration
  filter(
    brand %in% c("B1", "B2"),
    item  %in% c("1", "2")
  ) |>
  # Prepare the data for modeling
  par_model_mat(
    # Specify the model formula
    QTY ~ PROMO + brand*item,
    # Specify the time variable to ensure correct ordering
    time = DATE,
    # Specify grouping variables to handle hierarchical structure
    groups = c(brand, item),
    # Number of lags for the autoregressive component
    nlag = 2
  )

# Inspect the prepared model specification
model_spec
```

    PAR model specification:
     - QTY ~ lag1 + lag2 + (PROMO + brand * item) 
     - design dim: 7184 x 7 
     - across 4 group(s)
     - using 2 lag point(s) and 5 covariate(s)

``` r
# Step 2: Fit the model using MCMC
res <- fit_par_mcmc(
  model_spec$Y,
  model_spec$X,
  q = length(model_spec$beta)
)

# View the estimated parameters
res$beta  # Autoregressive coefficients
```

    [1] 0.2262182 0.1030291

``` r
res$gamma # Covariate effects
```

    [1]  1.387327  1.438974 -1.973517 -0.789424  1.823724

``` r
# Step 3: Fit the model using BFGS optimization
res2 <- fit_par_bfgs(model_spec)

# View the estimated parameters from BFGS
res2$beta  # Autoregressive coefficients
```

         lag1      lag2 
    0.2269197 0.1030823 

``` r
res2$gamma # Covariate effects
```

      (Intercept)         PROMO       brandB2         item2 brandB2:item2 
        1.3899075     1.4335929    -1.9683958    -0.7882043     1.8212733 

For examples of further analysis, please see either the vignette or
example analysis/scripts within the `inst/` directory of the github.

If `pastasales` is installed, you may view these files with the
following R code:

``` r
library(pastasales)
browseVignettes("pastasales")
# Note some additional packages may be required (listed in Suggests)
file.edit(system.file("analysis/bayes_analysis.R", package = "pastasales"))
file.edit(system.file("analysis/eda.R", package = "pastasales"))
file.edit(system.file("analysis/RandomForest.r", package = "pastasales"))
file.edit(system.file("analysis/TS_model.R", package = "pastasales"))
file.edit(system.file("examples/em-results.r", package = "pastasales"))
```
