#' Best Beta
#'
#'
#' @description Best parameter for penality
#'
#' @param data vector of data points
#'
#' @return a number
#' @export
#'
#' @examples
#' best_beta(c(0,1,1,2,2,3))
best_beta <- function(data){
  n = length(data)
  return (2*var(data)*log(n))
}
