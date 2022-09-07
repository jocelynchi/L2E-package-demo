#' L2E multivariate regression - PG
#'
#' \code{l2e_regression} performs L2E multivariate regression via block coordinate descent
#' with proximal gradient for updating both beta and tau.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param tau Initial precision estimate
#' @param b Initial vector of regression coefficients
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (vector) and tau (scalar),
#' the number of outer block descent iterations until convergence (scalar),
#' and the number of inner iterations per outer iteration for updating beta and tau (vectors)
#' @importFrom stats sd
#' @export
#' @examples
#' # Bank data example
#' y <- bank$y
#' X <- as.matrix(bank[,1:13])
#' X0 <- as.matrix(cbind(rep(1,length(y)), X))
#' tauinit <- 1/mad(y)
#' binit <- matrix(0, 14, 1)
#'
#' sol <- l2e_regression(y, X0, binit, tauinit)
#' r <- y - X0 %*% sol$beta
#' ix <- which(abs(r) > 3/sol$tau)
#' l2e_fit <- X0 %*% sol$beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
l2e_regression <- function(y,X,b,tau,max_iter=1e2,tol=1e-4,Show.Time=TRUE) {

  if (tau <= 0) stop("Entered non-positive initial tau")

  time <- proc.time()

  # save the inner loop iters
  iter_beta <- rep(0, max_iter)
  iter_tau <- rep(0, max_iter)

  n <- nrow(X); p <- ncol(X)
  if (p >= n) stop('Current implementation can only handle n > p')
  QRF <- qr(X)
  sd_y <- as.numeric(sd(y))

  for (i in 1:max_iter) {
    b_last <- b
    tau_last <- tau

    res_b <- update_beta_qr(y,X,QRF,tau_last,b_last)
    b <- res_b$beta
    iter_beta[i] <- res_b$iter

    r <- y - X%*%b
    res_tau <- update_tau_R(r,tau_last,sd_y)
    tau <- res_tau$tau
    iter_tau[i] <- res_tau$iter

    A <- norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))
    B <- abs(tau_last-tau) < tol*(1 + tau_last)
    if (A & B) break
  }
  if(Show.Time) print(proc.time() - time)

  return(list(beta=b,tau=tau,iter=i,iter_beta=iter_beta[1:i],iter_tau=iter_tau[1:i]))
}
