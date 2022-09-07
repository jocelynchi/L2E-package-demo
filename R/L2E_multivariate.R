#' L2E multivariate regression
#'
#' \code{L2E_multivariate} performs multivariate regression under the L2 criterion. Available methods include proximal gradient descent (PG) and majorization-minimization (MM).
#'
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param method Available methods include PG and MM. MM by default.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (vector) and tau (scalar),
#' the number of outer block descent iterations until convergence (scalar),
#' and the number of inner iterations per outer iteration for updating beta (vector) and tau or eta (vector)
#' @export
#' @examples
#' # Bank data example
#' y <- bank$y
#' X <- as.matrix(bank[,1:13])
#' X0 <- as.matrix(cbind(rep(1,length(y)), X))
#'
#' tau <- 1/mad(y)
#' b <- matrix(0, 14, 1)
#'
#' # MM method
#' sol_mm <- L2E_multivariate(y, X0, b, tau)
#' r_mm <- y - X0 %*% sol_mm$beta
#' ix_mm <- which(abs(r_mm) > 3/sol_mm$tau)
#' l2e_fit_mm <- X0 %*% sol_mm$beta
#'
#' # PG method
#' sol_pg <- L2E_multivariate(y, X0, b, tau, method="PG")
#' r_pg <- y - X0 %*% sol_pg$beta
#' ix_pg <- which(abs(r_pg) > 3/sol_pg$tau)
#' l2e_fit_pg <- X0 %*% sol_pg$beta
#'
#' plot(y, l2e_fit_mm, ylab='Predicted values', main='MM', pch=16, cex=0.8) # MM
#' points(y[ix_mm], l2e_fit_mm[ix_mm], pch=16, col='blue', cex=0.8) # MM
#' plot(y, l2e_fit_pg, ylab='Predicted values', main='PG', pch=16, cex=0.8) # PG
#' points(y[ix_pg], l2e_fit_pg[ix_pg], pch=16, col='blue', cex=0.8) # PG
#'
L2E_multivariate <- function(y,X,beta,tau,method="MM",max_iter=1e2,tol=1e-4,Show.Time=TRUE) {

  if (tau <= 0) stop("Entered non-positive initial tau")

  if (method == "MM") {
    l2e_regression_MM(y=y, X=X, beta=beta, tau=tau, max_iter=max_iter, tol=tol, Show.Time=Show.Time)
  }
  else if (method == "PG") {
    l2e_regression(y=y, X=X, b=beta, tau=tau, max_iter=max_iter, tol=tol, Show.Time=Show.Time)
  }
  else {
    stop("Method unavailable. Available methods are PG and MM.")
  }
}
