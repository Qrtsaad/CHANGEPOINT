#' Cost Gauss
#'
#'
#' @description Cost function for Gaussian models
#'
#' @param trajectoire vector of trajectoire from data points
#'
#' @return a number
#' @export
#'
#' @examples
#' cost_gauss(c(0,0,1,1,2,2,3,4,5))
cost_gauss <- function(trajectoire) {
  n = length(trajectoire)
  if (n == 1) {res = 0} else {res = (n-1)*var(trajectoire)}
  return (res)
}


#' Cost Poisson
#'
#' @description Cost function for Poisson models
#'
#' @param trajectoire vector of trajectoire from data points
#'
#' @return a number
#' @export
#'
#' @examples
#' cost_poiss(c(1,2,3))
#' cost_poiss(c(1,1,1,2,2,2))
cost_poiss <- function(trajectoire) {
  n <- length(trajectoire)
  moy <- mean(trajectoire)
  slogfact <- sum(log(factorial(trajectoire)))
  return (2*n*(moy - log(moy)*moy + slogfact))

}


