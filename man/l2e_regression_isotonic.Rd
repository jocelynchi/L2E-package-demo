% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l2e_regression_isotonic.R
\name{l2e_regression_isotonic}
\alias{l2e_regression_isotonic}
\title{LTE isotonic regression}
\usage{
l2e_regression_isotonic(y, b, tau, max_iter = 100, tol = 1e-04)
}
\arguments{
\item{y}{Response vector}

\item{b}{Initial vector of regression coefficients}

\item{tau}{Initial precision estimate}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Relative tolerance}
}
\value{
Returns a list object containing the estimates for beta (vector) and tau (scalar)
}
\description{
\code{l2e_regression_isotonic} Performs L2E isotonic regression via inexact coordinate descent.
}
\examples{

set.seed(12345)
n <- 200
tau <- 1
x <- seq(-2.5, 2.5, length.out=n)
f <- x^3
y <- f + (1/tau)*rnorm(n)

# Clean Data
plot(x, y, pch=16)
lines(x, f, col='blue', lwd=3)

tau <- 1
b <- y
sol <- l2e_regression_isotonic(y, b, tau)

plot(x, y, pch=16)
lines(x, f, col='blue', lwd=3)
iso <- gpava(1:n, y)$x
lines(x, iso, col='red', lwd=3)
lines(x, sol$beta, col='green', lwd=3)

# Contaminated Data
ix <- 0:9
y[45 + ix] <- 14 + rnorm(10)

plot(x, y, pch=16)
lines(x, f, col='blue', lwd=3)

tau <- 1
b <- y
sol <- l2e_regression_isotonic(y, b, tau)

plot(x, y, pch=16)
lines(x, f, col='blue', lwd=3)
iso <- gpava(1:n, y)$x
lines(x, iso, col='red', lwd=3)
lines(x, sol$beta, col='green', lwd=3)

}
