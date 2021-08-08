#' Function NegBin
#'
#' @param x variable of function
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#'
#' @return a value or vector
#' @export
#'
#' @examples
#' J(0.5, 20, 50, 1)
#' J(c(0.2,0.3,0.4,0.5), 20, 50, 1)
J <- function(x, A=20, B=30, C=1)
{
  return ((A/C)*log(x) + (B/C)*log(1-x) + 1)
}

#' Function derivative NegBin
#'
#' @param x variable of function
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#'
#' @return a value or vector
#' @export
#'
#' @examples
#' dJ(0.5, 20, 50, 1)
#' dJ(c(0.2,0.3,0.4,0.5), 20, 50, 1)
dJ <- function(x,A=20,B=30,C=1)
{
  return (A/(C*x) - B/(C*(1-x)))
}


#' good_values
#'
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#'
#' @return a boolean
#' @export
#'
#' @examples
#' good_values(-10,5,100)
#' good_values(10,39,24)
good_values <- function(A,B,C)
{
  res <- TRUE
  if (A < 0 | B < 0 | C < 0){res <- FALSE}
  return (res)
}


#' exist_zero
#'
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#'
#' @return a boolean
#' @export
#'
#' @examples
#' exist_zero(-10,5,100)
#' exist_zero(10,39,24)
exist_zero <- function(A,B,C)
{
  res <- TRUE
  if (J(A/(A+B),A,B,C) < 0){res <- FALSE}
  return (res)
}


#' xABC
#'
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#'
#' @return 3 values of initial guess for root-finding algorithms
#' @export
#'
#' @examples
#' xABC(100,1,100)
#' xABC(500,10,200)
#' xABC(2000,3,1500)
xABC <- function(A,B,C)
{

  if (A >= 2*C)
  {
    x0 <- runif(1,min = 0.7, max = 0.9)
    x1 <- runif(1,min = 0.8, max = 1)
    x2 <- (x0+x1)/2
  }
  else
  {
    if (C <= 4*B)
    {
      x0 <- runif(1,min = 0.2, max = 0.4)
      x1 <- runif(1,min = 0.8, max = 1)
      x2 <- (x0+x1)/2
    }
    else
    {
      x0 <- runif(1,min = 0.6, max = 0.8)
      x1 <- runif(1,min = 0.7, max = 0.9)
      x2 <- (x0+x1)/2
    }
  }

  return(list(x0 = x0, x1 = x1, x2 = x2))
}


#' searchzero
#'
#' @param A parameter of function
#' @param B parameter of function
#' @param C parameter of function
#' @param tol desired precision
#' @param method name of the root-finding algorithm used
#' @param x0 initial guess for root-finding algorithms
#' @param x1 initial guess for root-finding algorithms
#' @param x2 initial guess for root-finding algorithms
#'
#' @return a list of values with the approximate solution, the number of iterations and the complexity of the method used
#' @export
#'
#' @examples
#' searchzero(A = 1000, B = 10, C = 100, method = "newton", x0 = 0.2)
#' searchzero(A = 500, B = 5, C = 1000, method = "secante")
searchzero <- function(A = 20, B = 30, C = 1, tol = 1e-10, method, x0 = xABC(A,B,C)$x0, x1 = xABC(A,B,C)$x1, x2 = xABC(A,B,C)$x2)
{
  allowed.method <- c("newton", "muller", "secante")
  if(!method %in% allowed.method){stop('type must be one of: ', paste(allowed.method, collapse=", "))}

  if(good_values(A,B,C) == FALSE){stop('the given values are incorrect ')}

  else if(exist_zero(A,B,C) == FALSE){stop('there is no solution with these values')}

  else
  {

    if (method == "newton")
    {
      x <- x0
      it <- 0
      tmp <- 0
      while (abs(J(x,A,B,C)) > tol)
      {
        tmp <- x
        x <- x - J(x,A,B,C)/dJ(x,A,B,C)
        it <- it+1
      }
      cp <- abs(x-tmp)
    }

    else if (method == "secante")
    {
      xm1 <- x0
      x <- x1
      it <- 0
      tmp <- 0
      while (abs(J(x,A,B,C)) > tol)
      {
        tmp <- x
        x <- x - J(x,A,B,C)*(x - xm1)/(J(x,A,B,C) - J(xm1,A,B,C))
        xm1 <- tmp
        it <- it+1
      }
      cp <- abs(x-tmp)
    }

    else if (method == "muller")
    {
      J <- function(x) (A/C)*log(x) + (B/C)*log(1-x) + 1
      dJ <- function(x) A/(C*x) - B/(C*(1-x))

      x <- muller(J, p0 = x0, p1 = x1, p2 = x2, tol = tol)$root
      it <- muller(J, p0 = x0, p1 = x1, p2 = x2, tol = tol)$niter
      cp <- muller(J, p0 = x0, p1 = x1, p2 = x2, tol = tol)$reltol
    }
  }

  return (list(x_star = x, iterations = it, complexity = cp))
}
