#' Beta update in L2E isotonic regression - PG
#'
#' \code{update_beta_isotonic} updates beta for L2E isotonic regression using PG
#'
#' @param y Response vector
#' @param b Current estimate for beta
#' @param tau Current estimate for tau
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#' @importFrom isotone gpava
#'
update_beta_isotonic <- function(y,b,tau,max_iter=1e2,tol=1e-4) {
  n <- length(y)
  for (i in 1:max_iter) {
    b_last <- b
    r <- y - b
    w <- exp(-0.5* (tau*r)**2 )
    z <- w*y + (1-w)*b
    b <- gpava(1:n,z)$x
    if ( norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))) break
  }
  return(list(beta=b,iter=i))
}
