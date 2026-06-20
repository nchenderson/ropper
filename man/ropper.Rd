\name{ropper}
\alias{ropper}
\title{
Ranking-targeted estimation of covariate-adjusted rankings
}
\description{
Computes ranking of cluster-specific estimates after performing a covariate
adjustment.
}
\usage{
ropper(Y, X, ses, tau.sq = c("REML", "kNN"), H=1, 
       opt.method=c("MM", "optim"), trim.quant=0,
       control = list())
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
  \item{tau.sq}{
Method for estimating the random effects variance. This can be either the 
restricted maximum likelihood estimation (REML) \code{REML} or k-nearest 
neighbors \code{kNN}.
}
  \item{H}{
The degree of the approximation used for the unbiased risk estimator. This must 
be an integer greater than or equal to 1. 
}
  \item{opt.method}{
The optimization method used for computing the regression coefficient estimates. 
This can be either an MM algorithm or the use of the \code{optim} function.       
}
\item{trim.quant}{
The proportion of extreme estimates to exclude when estimating the random effects
variance. Cluster-specific estimates lower than the \code{trim.quant} quantile will be 
set equal to the \code{trim.quant} quantile, and estimates larger than the
1 - \code{trim.quant} quantile will be set equal to the 1 - \code{trim.quant} quantile.
This must be a number greater than or equal to 0 and less than 0.5.
}
\item{control}{A list of control parameters specifying any changes to
      default values of the optimization procedures. Full names of control 
      list elements must be specified, otherwise, user-specifications are ignored.
	    See *Details*.}
}
\details{
Default values of \code{control} are:
    \code{maxiter=500},
    \code{tol=1e-06}

\describe{
	\item{\code{maxiter}}{An integer denoting the maximum limit on the number
	  of MM iterations or \code{optim} iterations.  Default value is 500.}
  \item{\code{tol}}{A small, positive scalar that determines when iterations
	  should be terminated. Iterations are terminated when
	  \eqn{||\beta_{k+1} - \beta_k|| < \textrm{tol}}{|| \beta[k+1] - \beta[k] || < tol}.
	  Default is \code{1.e-06}.}
}
}
\value{
A list with the following components:
\item{coefficients}{The length-p vector of ranking-targeted regression coefficient estimates}
\item{ppep}{The length-K vector of estiamted population posterior expected percentiles}
\item{objfn.track}{A vector containing the value of the objective function at each iteration}
}
\references{
Henderson, N. C., & Hartman, N. (2025). Targeted Parameter Estimation for Robust Empirical Bayes Ranking. 
\emph{arXiv preprint}. \url{https://arxiv.org/abs/2511.16530}
}
\author{
Nicholas C. Henderson and Nicholas Hartman
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#######################################################
###  Example of using ropper with simulated data  #####
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
stderrs_sq <- 1/cluster.sizes ## define standard errors squared

rfit <- ropper(Y=yy, X=XX, ses=stderrs_sq)

## Regression coefficient estimates
rfit$coefficients

## ROPPER percentile estimates for first 5 units
rfit$ppep[1:5]

## Plot true random effect rankings vs. ROPPer rankings:
true_ranks <- rank(theta.vec)
ropper_ranks <- rank(rfit$ppep)

plot(true_ranks, ropper_ranks, las=1,
     xlab="True RE Ranking", ylab="ROPPER Ranking")
abline(0, 1)
}
\keyword{ regression }
\keyword{ robust }
\concept{ empirical Bayes }
\concept{ loss-function driven estimation }
\concept{ unbiased risk estimation }
\concept{ ranking }
\concept{ mixed models }
