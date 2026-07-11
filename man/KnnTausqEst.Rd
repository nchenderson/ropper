\name{KnnTausqEst}
\alias{KnnTausqEst}
\title{
Nearest neighbors estimation of random effects variance
}
\description{
Computes an estimate of the random effects variance using a nonparametric, 
k nearest neighbors approach. In this approach, the standard
errors for the cluster-level effect estimates are assumed to be fixed.
}
\usage{
KnnTausqEst(Y, X, ses, trim.quant = 0, k = 1)
}
\arguments{
  \item{Y}{
The K x 1 vector of cluster-specific effect estimates.
}
  \item{X}{
The K x p design matrix containing cluster-specific covariates.
}
  \item{ses}{
The K x 1 vector sampling standard errors squared.
}
  \item{trim.quant}{
The proportion of extreme estimates to exclude when estimating the random effects
variance. Cluster-specific estimates lower than the \code{trim.quant} quantile will be 
set equal to the \code{trim.quant} quantile, and estimates larger than the
1 - \code{trim.quant} quantile will be set equal to the 1 - \code{trim.quant} quantile.
This must be a number greater than or equal to 0 and less than 0.5.
}
\item{k}{
The number of neighbors to include in the k nearest neighbors procedure. This should
be an integer both greater than or equal to 1 and less than or equal to K.
}
}
\details{
Calculation of the nearest neighbor estimator of the random effects variance involves a random
split of the data into a training set \eqn{D_{1}} and a test set \eqn{D_{2}} of approximately equal sizes.
The one-nearest neighbor estimator (i.e., when \code{k=1}) of \eqn{\tau^{2}} is 
then calculated as
\eqn{\hat{\tau}^{2}_{1-NN} = \frac{1}{K}\sum_{j=1}^{K} Y_{j}^{2} - \frac{1}{K_{T}}\sum_{j \in D_{2}} Y_{j}Y_{j}' 
- \frac{1}{K_{T}}\sum_{j \in D_{2}} ses_{j}^{2}}, where \eqn{K_{T}} is the number of observations in the test
set \eqn{D_{2}}, and \eqn{Y_{j}'} is the effect estimate associated 
with the nearest neighbor \eqn{\mathbf{x}_{j}'} 
of \eqn{\mathbf{x}_{j}} in the training set.
}
\value{
The estimate of the random effects variance. A vector of length 1.
}
\references{
Devroye, L., Györfi, L., Lugosi, G., & Walk, H. (2018). A nearest neighbor 
estimate of the residual variance. Electronic Journal of Statistics, 12(1), 
1752–1778. \doi{10.1214/18-EJS1438}
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
See also \code{\link{REMLTausqEst}}, \code{\link{ropper}}
}
\examples{
##########################################################
###  Example of kNN estimation with simulated data  #####
##########################################################
n <- 500   ## 500 units
## True tau.sq is 0.5
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

## Run Knn estimation with a 1 nearest neighbor approach
tausq_hat_knn <- KnnTausqEst(Y=yy, X=XX, ses=stderrs_sq, k=1)
tausq_hat_knn
}
\keyword{ regression }
\keyword{ robust }
