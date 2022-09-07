#' Beta update in L2E convex regression - MM
#'
#' \code{update_beta_MM_convex} updates beta for L2E convex regression using MM
#'
#' @param y Response
#' @param beta Initial vector of regression coefficients
#' @param tau Precision estimate
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
#' @importFrom cobs conreg
#'

update_beta_MM_convex <- function(y,beta,tau,max_iter=1e2,tol=1e-4) {

  n <- length(y)
  for (i in 1:max_iter) {
    beta_last <- beta
    # Compute the weights
    r <- y - beta_last
    w <- exp(-0.5* (tau*r)**2 )

    # Now solve a weighted convex regression,
    # cobs package has a function conreg, solving \sum_{i=1}^n w_ii (y_i - \beta_i)^2
    beta <- conreg(x=1:n, y, sqrt(w), convex = TRUE, method = "Duembgen06_R")$yf
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  return(list(beta=beta,iter=i))
}
