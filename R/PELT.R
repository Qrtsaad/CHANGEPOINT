
#' Pruned Exact Linear Time
#'
#'
#' @description Pruned Exact Linear Time Algorithm
#'
#' @param data vector of data points
#' @param cost a number
#' @param beta a number
#'
#' @return a vector of changepoints, global cost
#' @export
#'
#' @examples
#' myPELT(c(0,0,1,1,0,0,0), beta = 0.00001)
#' myPELT(c(rnorm(50, mean = 0, sd = 0), rnorm(50, mean = 10, sd = 0)), beta = 1)
#' myPELT(data_generator(25, chpts = c(10,20), means = c(20,0,20), type = "gauss"), beta = 5)

myPELT <- function(data, cost = "gauss", beta = best_beta(data))
{
  allowed.cost <- c("gauss", "poisson")
  if(!cost %in% allowed.cost){stop('type must be one of: ', paste(allowed.cost, collapse=", "))}

  if (cost == "gauss") {cost_f <- cost_gauss}
  else if (cost == "poisson") {cost_f <- cost_poiss}

  n <- length(data)

  cp <- rep(0,n)
  Q <- rep(0,n)
  R <- NULL

  for (t in 2:n)
  {
    val_min <- cost_f(data[1:t])
    arg_min <- 0
    R <- c(R,t)

    for (s in R)
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

    for (s in R){
      #     print(s)
      #     print(t)
      #     print(Q[s] + cost_f(data[(s+1):t]))
      if (Q[t] <= Q[s-1] + cost_f(data[s:t])){
        R <- R[R!= s]}}
  }

#backtracking
  v <- cp[n]
  P <- cp[n]

  while (v > 0)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(list(changepoints = P, mean = NULL, globalCost = Q[n] - length(P)*beta))
}

