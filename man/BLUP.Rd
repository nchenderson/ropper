\name{BLUP}
\alias{BLUP}
\title{
Posterior expected means
}
\description{
Computes estimates of cluster-level posterior expected means given an estimate of 
the regression coefficients and an estimate of the random effects variance. For
a fixed value of the random effects variance, the collection of posterior expected
mean estimates is commonly referred to as the "best linear unbiased predictor",
or BLUP.
}
\usage{
BLUP(beta.hat, Y, X, ses, tau.sq)
}
\arguments{
  \item{beta.hat}{
The 1 x p vector of regression coefficient estimates.
}
  \item{Y}{
The K x 1 vector of cluster-specific effect estimates.
}
  \item{X}{
The K x p design matrix containing cluster-specific covariates.
}
  \item{ses}{
The K x 1 vector sampling standard errors squared.
}
  \item{tau.sq}{
Estimate of the random effects variance. This must 
be a scalar greater than 0. 
}
}
\details{
The BLUP for cluster k is defined as
\eqn{E(v_{k} | Y_{k} ) = B_{k}( Y_{k} - X_{k}^{T}\hat{\beta}) }, where
\eqn{\hat{\beta}} is the vector of regression coefficient estimates and \eqn{B_{k}}
is defined as \eqn{B_{k} = \tau^{2}/(\tau^{2} + ses_{k}) }
}
\value{
A K x 1 vector of posterior mean estimates.
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
See also \code{\link{PEPP}}.
}
\examples{
#######################################################
###  Example of using PEPP with simulated data  #####
#######################################################
n <- 50
tau.sq <- 0.5
cluster.sizes <- 1 + rpois(n, lambda=5)

xx <- rnorm(n) 
theta.vec <- rnorm(n, sd=sqrt(tau.sq))
## True fixed-effect function is pnorm(x)
mu.vec <- pnorm(xx) + theta.vec 
yy <- mu.vec + rnorm(n)/sqrt(cluster.sizes)

## Fit model where it is assumed that fixed-effect
## function is beta0 + beta1*x
XX <- cbind(rep(1, n), xx)
## define standard errors squared
stderrs_sq <- 1/cluster.sizes 


## Estimate regression coefficients by maximum likelihood 
## estimation assuming that tau.sq = 0.5
beta.mle <- solve(crossprod(XX, (stderrs_sq + 0.5)*XX), 
                  crossprod(XX, (stderrs_sq + 0.5)*yy)) 
## Extract the posterior means associated with the maximum 
## likelihood estimate:
blup.mle <- BLUP(beta.mle, Y=yy, X=XX, ses=stderrs_sq, 
                 tau.sq=0.5)

## Estimate regression coefficients using ROPPER assuming
## that tau.sq = 0.5
rfit <- ropper(Y=yy, X=XX, ses=stderrs_sq, tau.sq=0.5)

## Extract posterior means associated with the ROPPER estimate
blup.rop <- BLUP(rfit$coefficients, Y=yy, X=XX, 
                 ses=stderrs_sq, tau.sq=0.5)

plot(blup.mle, blup.rop, 
     xlab="MLE Posterior Means", 
     ylab="ROPPER Posterior Means")
}
\keyword{ regression }
\concept{ empirical Bayes }
\concept{ shrinkage }
