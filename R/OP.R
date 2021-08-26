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
#' myOP(data_generator(25, chpts = c(10,20), means = c(20,0,20), type = "gauss"), beta = 5)
myOP <- function(data, cost = "gauss", beta = best_beta(data), eps = 1e-6)
{
  allowed.cost <- c("gauss", "poisson", "negbin")
  if(!cost %in% allowed.cost){stop('type must be one of: ', paste(allowed.cost, collapse=", "))}

  if (cost == "gauss") {cost_f <- cost_gauss}
  else if (cost == "poisson") {cost_f <- cost_poiss}
  else if (cost == "negbin")
  {
    phi <- mean(data)^2/(sd(data) - mean(data))
    data <- data/phi
    data[data==0] <- eps/(1-eps)
    cost_f <- cost_negbin
  }

  n <- length(data)

  cp <- rep(0,n)
  Q <- rep(0,n)

  for (t in 2:n)
  {
    val_min <- cost_f(data[1:t])
    arg_min <- 0
    for (s in 2:t)
    {
      a <- Q[s-1] + cost_f(data[s:t]) + beta
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

  return(list(changepoints = P, globalCost = Q[n] - length(P)*beta))
}

