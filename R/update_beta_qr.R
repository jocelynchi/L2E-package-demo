#' Beta update in L2E multivariate regression - PG
#'
#' \code{update_beta_qr} updates beta for L2E multivariate regression via a QR solve
#'
#' @param y Response vector
#' @param X Design matrix
#' @param QRF QR factorization object for X (obtained via `QRF=qr(X)`)
#' @param b Current estimate for beta
#' @param tau Current estimate for tau
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#'
update_beta_qr <- function(y,X,QRF,tau,b,max_iter=1e2,tol=1e-4) {
  tau <- as.numeric(tau)
  for (i in 1:max_iter) {
    b_last <- b
    Xb <- X %*% b
    r <- y - Xb
    w <- exp(-0.5* (tau*r)**2 )
    u <- w*r + Xb
    b <- qr.solve(QRF, u)
    if (norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))) break
  }
  return(list(beta=b,iter=i))
}
