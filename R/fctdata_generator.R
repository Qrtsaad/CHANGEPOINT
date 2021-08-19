#' Make changepoints
#'
#' @param n an integer
#'
#' @return a vector of changepoints for data_generator function
#' @export
#'
#' @examples
make_chpts <- function(n)
{
  if (n<0){stop('n must be non-negative')}
  else if (n <= 100){res <- seq(from = 0, to = n-1, by = n/2)}
  else if (n <= 200){res <- seq(from = 0, to = n-1, by = n/4)}
  else if (n <= 500){res <- seq(from = 0, to = n-1, by = n/5)}
  else if (n <= 1000){res <- seq(from = 0, to = n-1, by = n/10)}
  else {res <- seq(from = 0, to = n-1, by = n/20)}
  return (res)
}

#' Make means
#'
#' @param n an integer
#'
#'
#' @return a vector of means for data_generator function
#' @export
#'
#' @examples
make_means <- function(n)
{
  res <- NULL
  tmp <- 0
  for (i in 1:(n+1)){
    rand <- sample(0:10, 1)
    while (rand == tmp){rand <- sample(1:10, 1)}
    tmp <- rand
    res <- c(res,rand)
  }
  return (res)
}
