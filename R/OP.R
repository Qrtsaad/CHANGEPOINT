#' Optimal Partitioning
#'
#'
#' @description Optimal Partitioning Algorithm
#'
#' @param data vector of data points
#' @param cost type of my cost (gauss or poisson)
#' @param beta penalty value
#'
#' @return a vector of changepoints, a number for the complexity (cost of computations)
#' @export
#'
#' @examples
#' myOP(c(0,0,1,1,0,0,0), beta = 0.00001)
#' myOP(c(rnorm(50, mean = 0, sd = 0), rnorm(50, mean = 10, sd = 0)))
myOP <- function(data, cost = "gauss", beta = best_beta(data))
{
  allowed.cost <- c("gauss", "poisson")
  if(!cost %in% allowed.cost){stop('type must be one of: ', paste(allowed.cost, collapse=", "))}

  if (cost == "gauss") {cost_f <- cost_gauss}
  else if (cost == "poisson") {cost_f <- cost_poiss}

  n <- length(data)

  cp <- rep(0,n)
  Q <- rep(0,n)

  for (t in 2:n)
  {
    val_min <- cost(data[1:t])
    arg_min <- 0
    for (s in 2:t)
    {
      a <- Q[s-1] + cost(data[s:t]) + beta
      if (a < val_min)
      {
        val_min <- a
        arg_min <- s - 1
      }
    }
    Q[t] <- val_min
    cp[t] <- arg_min
  }

# backtracking
  v <- cp[n]
  P <- cp[n]

  while (v > 0)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(list(changepoints = P, means = NULL, globalCost = Q[n] - length(P)*beta))
}

