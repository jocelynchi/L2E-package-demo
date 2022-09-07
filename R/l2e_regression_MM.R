#' L2E multivariate regression - MM
#'
#' \code{l2e_regression_MM} performs L2E multivariate regression via block coordinate descent
#' with MM for updating beta and modified Newton for updating tau.
#'
#' @param y response
#' @param X Design matrix
#' @param beta initial vector of regression coefficients
#' @param tau initial precision estimate
#' @param max_iter maximum number of iterations
#' @param tol relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (vector) and tau (scalar),
#' the number of outer block descent iterations until convergence (scalar),
#' and the number of inner iterations per outer iteration for updating beta and eta (vectors)
#' @export
#' @examples
#' # Bank data example
#' y <- bank$y
#' X <- as.matrix(bank[,1:13])
#' X0 <- as.matrix(cbind(rep(1,length(y)), X))
#' tau <- 1/mad(y)
#' b <- matrix(0, 14, 1)
#'
#' sol <- l2e_regression_MM(y, X0, b, tau)
#' r <- y - X0 %*% sol$beta
#' ix <- which(abs(r) > 3/sol$tau)
#' l2e_fit <- X0 %*% sol$beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
l2e_regression_MM <- function(y,X,beta,tau,max_iter=1e2,tol=1e-4,Show.Time=TRUE) {

  if (tau <= 0) stop("Entered non-positive initial tau")

  time <- proc.time()

  # save the inner loop iters
  iter_beta <- rep(0, max_iter)
  iter_eta <- rep(0, max_iter)

  for (i in 1:max_iter) {
    # update beta
    beta_last <- beta
    res_beta <- update_beta_MM_ls(y,X,beta_last,tau,max_iter=1e2,tol=1e-4)
    beta <- res_beta$beta
    iter_beta[i] <- res_beta$iter

    # update tau/eta
    r <- y - X%*%beta
    eta_last <- log(tau)  # get eta as in line 9
    res_eta <- update_eta_bktk(r, eta_last, tol=tol) # update eta as in line 10-12
    eta <- res_eta$eta
    tau <- exp(eta) # update tau as in line 13
    iter_eta[i] <- res_eta$iter


    A <- norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))
    B <- abs(eta_last-eta) < tol*(1 + abs(eta_last))  # tau_last could be negative, so need abs() here
    if (A & B) break
  }
  if(Show.Time) print(proc.time() - time)

  return(list(beta=as.vector(beta),tau=tau, iter=i,
              iter_beta = iter_beta[1:i], iter_eta = iter_eta[1:i]))
}
