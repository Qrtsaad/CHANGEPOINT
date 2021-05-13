#' Functional Pruned Optimal Partitioning 1D
#'
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm
#'
#' @param data vector of data points
#' @param cost a number
#' @param beta a number
#'
#' @return a vector of changepoints, global cost
#' @export
#'
#' @examples
#' myFPOP1D(c(0,0,1,1,0,0,0), beta = 0.00001)
#' myFPOP1D(c(rnorm(50, mean = 0, sd = 0), rnorm(50, mean = 10, sd = 0)), beta = 1)
myFPOP1D <- function(data, cost = "gauss", beta = best_beta(data))
{
  allowed.cost <- c("gauss", "poisson")
  if(!cost %in% allowed.cost){stop('type must be one of: ', paste(allowed.cost, collapse=", "))}

  if (cost == "gauss") {cost_f <- cost_gauss}
  else if (cost == "poisson") {cost_f <- cost_poiss}

  n <- length(data)
  m <- -beta

  tau <- rep(0,n)
  v <- NULL

  for (t in 2:n)
  {
    getV <- getValue(v, t-1, beta)
    u <- list(label = t, value = getV + beta, position = data[t])
    v[t-1] <- list(v[t-1], u)

    v[t] <- v[t-1] # adding the new data point yt
    thetastar <- mean(data[1:t])
    v[t] <- list(v[t], quad(data, t, t, beta, data[t])) # adding intersection point new quadratic qtt vs old ones qti


    v[t] <- sort(v[1:t]) # sort with R function
    v[t] <- prune(v[1:t]) # pruning

    tau[t] <- getLabel(data, t, beta) # label of the smallest element

    # backtracking
    s <- tau[n]
    P <- tau[n]

    while (s)
    {
      P <- c(P, cp[s])
      s <- cp[s]
    }
    P <- rev(P)[-1]

  }


  return(list(changepoints = P, mean = NULL, globalCost = v[n] - length(P)*beta))
}
