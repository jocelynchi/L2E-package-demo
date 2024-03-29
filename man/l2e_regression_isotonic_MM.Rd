% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l2e_regression_isotonic_MM.R
\name{l2e_regression_isotonic_MM}
\alias{l2e_regression_isotonic_MM}
\title{L2E isotonic regression - MM}
\usage{
l2e_regression_isotonic_MM(
  y,
  beta,
  tau,
  max_iter = 100,
  tol = 1e-04,
  Show.Time = TRUE
)
}
\arguments{
\item{y}{Response vector}

\item{beta}{Initial vector of regression coefficients}

\item{tau}{Initial precision estimate}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Relative tolerance}

\item{Show.Time}{Report the computing time}
}
\value{
Returns a list object containing the estimates for beta (vector) and tau (scalar),
the number of outer block descent iterations until convergence (scalar),
and the number of inner iterations per outer iteration for updating beta and eta (vectors)
}
\description{
\code{l2e_regression_isotonic_MM} performs L2E isotonic regression via block coordinate descent
with MM for updating beta and modified Newton for updating tau.
}
\examples{
set.seed(12345)
n <- 200
tau <- 1
x <- seq(-2.5, 2.5, length.out=n)
f <- x^3
y <- f + (1/tau)*rnorm(n)

## Clean
plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
lines(x, f, lwd=3)

tau <- 1
b <- y
sol <- l2e_regression_isotonic_MM(y, b, tau)

plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
lines(x, f, lwd=3)
iso <- isotone::gpava(1:n, y)$x
lines(x, iso, col='blue', lwd=3)
lines(x,sol$beta,col='red',lwd=3)

## Contaminated
ix <- 0:9
y[45 + ix] <- 14 + rnorm(10)

plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
lines(x, f, lwd=3)

tau <- 1
b <- y
sol <- l2e_regression_isotonic_MM(y,b,tau)

plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
lines(x, f, lwd=3)
iso <- isotone::gpava(1:n, y)$x
lines(x, iso, col='blue', lwd=3)
lines(x,sol$beta,col='red',lwd=3)

}
