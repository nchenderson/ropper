\name{PEPP}
\alias{PEPP}
\title{
Posterior expected population percentiles
}
\description{
Computes estimates of cluster-level posterior expected population percentiles given an estimate of 
the regression coefficients and an estimate of the random effects variance.
}
\usage{
PEPP(beta.hat, Y, X, ses, tau.sq)
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
The posterior expected population percentile for cluster k is defined as
\eqn{E(\rho_{k} | Y_{k} ) = \Phi[ V_{k}( Y_{k} - X_{k}^{T}\hat{\beta})] }, where
\eqn{\Phi} is the cumulative distribution of a standard normal distribution, 
\eqn{\hat{\beta}} is the vector of regression coefficient estimates, and \eqn{V_{k}}
is defined as \eqn{V_{k} = \tau^{2}/[(\tau^{2} + ses_{k})\sqrt{2ses_{k} + \tau^{2}}] }.
The parameter \eqn{\rho_{k}} is the "population percentile" which is defined
as \eqn{\rho_{k} = \Phi(v_{k}/\tau)}, where \eqn{v_{k}} is the random effect
for the cluster-k effect estimate \eqn{Y_{k}}.
}
\value{
A K x 1 vector of posterior expected population percentile estimates.
}
\references{
Henderson, N. C., & Hartman, N. (2025). Targeted Parameter Estimation for Robust Empirical Bayes Ranking. 
\emph{arXiv preprint}. \url{https://arxiv.org/abs/2511.16530}
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
See also \code{\link{ropper}}, \code{\link{BLUP}}.
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
## Extract the PEPP associated with the maximum 
## likelihood estimate:
pepp.mle <- PEPP(beta.mle, Y=yy, X=XX, ses=stderrs_sq, 
                 tau.sq=0.5)

## Get PEPP using ROPPER
rfit <- ropper(Y=yy, X=XX, ses=stderrs_sq)
pepp_rop1 <- rfit$pepp
## This should match using PPEP with ROPPER estimates of 
## regression coefficients and tau.sq
pepp_rop2 <- PEPP(rfit$coefficients, Y=yy, X=XX, 
                  ses=stderrs_sq, tau.sq=rfit$tau.sq.hat)

all.equal(ppep_rop1, ppep_rop2)
}
\keyword{ regression }
\keyword{ robust }
\concept{ empirical Bayes }
\concept{ ranking }
