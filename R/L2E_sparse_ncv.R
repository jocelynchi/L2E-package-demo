#' Solution path of L2E sparse regression with existing penalization methods
#'
#' \code{L2E_sparse_ncv} computes the solution path of robust sparse regression under the L2 criterion. Available penalties include lasso, MCP and SCAD.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param b Initial vector of regression coefficients, can be omitted
#' @param tau Initial precision estimate, can be omitted
#' @param lambdaSeq A decreasing sequence of values for the tuning parameter lambda, can be omitted
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @return Returns a list object containing the estimates for beta (matrix) and
#' tau (vector) for each value of the tuning parameter lambda,
#' the run time (vector) for each lambda,
#' and the sequence of lambda used in the regression (vector)
#' @importFrom stats mad
#' @export
#' @examples
#' set.seed(12345)
#' n <- 100
#' tau <- 1
#' f <- matrix(c(rep(2,5), rep(0,45)), ncol = 1)
#' X <- X0 <- matrix(rnorm(n*50), nrow = n)
#' y <- y0 <- X0 %*% f + (1/tau)*rnorm(n)
#'
#' ## Clean Data
#' lambda <- 10^(-1)
#' sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' r <- y - X %*% sol$Beta
#' ix <- which(abs(r) > 3/sol$Tau)
#' l2e_fit <- X %*% sol$Beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
#' ## Contaminated Data
#' i <- 1:5
#' y[i] <- 2 + y0[i]
#' X[i,] <- 2 + X0[i,]
#'
#' sol <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' r <- y - X %*% sol$Beta
#' ix <- which(abs(r) > 3/sol$Tau)
#' l2e_fit <- X %*% sol$Beta
#'
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#'
L2E_sparse_ncv <- function(y,X,b,tau,lambdaSeq,penalty="MCP",max_iter=1e2,tol=1e-4,Show.Time=TRUE) {

  if(missing(b)){
    b <- double(ncol(X))  # initial beta
  }

  if(missing(tau)){
    tau <- 1/mad(y)   # initial tau
  }

  if (tau <= 0) stop("Entered non-positive initial tau")

  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of rho
  }


  Nlambda <- length(lambdaSeq)

  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)

  time <- proc.time()
  for (j in 1:Nlambda) {

    res <- l2e_regression_sparse_ncv(y, X, b, tau, lambda=lambdaSeq[j], penalty=penalty,
                                 max_iter=max_iter, tol=tol, Show.Time = FALSE)
    Beta[, j] <- b <- res$beta
    Tau[j] <- tau <- res$tau

  }
  runtime <-  proc.time() - time
  if(Show.Time) print(runtime)

  return(list(Beta=Beta, Tau=Tau, runtime=runtime, lambdaSeq=lambdaSeq))

}
