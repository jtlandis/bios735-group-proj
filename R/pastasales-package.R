#' BIOS735 Group Project (group 2)
#'
#' @description This package contains functions for analyzing pasta sales data. Included are a few rstan models for comparison.
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.3. https://mc-stan.org
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib pastasales, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstantools rstan_config
#' @importFrom rstan sampling
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr arrange
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr group_split
#' @importFrom dplyr if_all
#' @importFrom dplyr lag
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr starts_with
#' @importFrom dplyr ungroup
#' @importFrom magrittr %>%
#' @importFrom stats optim
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats setNames
#' @importFrom tibble tibble
## usethis namespace: end
NULL
