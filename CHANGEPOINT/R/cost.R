#' Cost Gauss
#'
#'
#' @description Cost function for Gaussian models
#'
#' @param v vector of data points
#'
#' @return the optimal cost value
#' @export
#'
#' @examples
#' cost_gauss(c(0,0,1,1))
cost_gauss <- function(v)
{
  n <- length(v)
  if (n == 1) {res <- 0} else {res <- (n-1)*var(v)}
  return (res)
}


#' Cost Poisson
#'
#' @description Cost function for Poisson models
#'
#' @param v vector of data points
#'
#' @return the optimal cost value
#' @export
#'
#' @examples
#' cost_poiss(c(0,0,0))
#' cost_poiss(c(1,1,1))
#' cost_poiss(c(1,1,1,2,2,2))
cost_poiss <- function(v)
{
  n <- length(v)
  moy <- mean(v)
  slogfact <- sum(log(factorial(v)))
  if (moy == 0){res <- 0} else {res <- n*moy*(1 - log(moy)) + slogfact}
  return (res)
}


