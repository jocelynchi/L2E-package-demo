#' LTE isotonic regression
#' 
#' \code{l2e_regression_isotonic} Performs L2E isotonic regression via inexact coordinate descent.
#' 
#' @param y Response vector
#' @param b Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the estimates for beta (vector) and tau (scalar)
#' @importFrom stats sd
#' @export
#' @examples
#' 
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' x <- seq(-2.5, 2.5, length.out=n)
#' f <- x^3
#' y <- f + (1/tau)*rnorm(n)
#' 
#' # Clean Data
#' plot(x, y, pch=16)
#' lines(x, f, col='blue', lwd=3)
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic(y, b, tau)
#' 
#' plot(x, y, pch=16)
#' lines(x, f, col='blue', lwd=3)
#' iso <- gpava(1:n, y)$x
#' lines(x, iso, col='red', lwd=3)
#' lines(x, sol$beta, col='green', lwd=3)
#' 
#' # Contaminated Data
#' ix <- 0:9
#' y[45 + ix] <- 14 + rnorm(10)
#' 
#' plot(x, y, pch=16)
#' lines(x, f, col='blue', lwd=3)
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic(y, b, tau)
#' 
#' plot(x, y, pch=16)
#' lines(x, f, col='blue', lwd=3)
#' iso <- gpava(1:n, y)$x
#' lines(x, iso, col='red', lwd=3)
#' lines(x, sol$beta, col='green', lwd=3)
#' 
l2e_regression_isotonic <- function(y,b,tau,max_iter=1e2,tol=1e-4) {
  sd_y <- sd(y)
  for (i in 1:max_iter) {
    b_last <- b
    tau_last <- tau
    b <- update_beta_isotonic(y,b,tau,max_iter=10,tol=1e-4)$beta 
    r <- y - b
    tau <- update_tau_R(r,tau_last, sd_y)$tau
    A <- norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))
    B <- abs(tau_last-tau) < tol*(1 + tau_last)
    if (A & B) break
  }
  return(list(beta=b,tau=tau))
}