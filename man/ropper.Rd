\name{ropper}
\alias{ropper}
\title{
Ranking-Optimized Population Posterior Expected Percentiles
}
\description{
Computes ranking of cluster-specific estimates after performing a covariate
adjustment.
}
\usage{
ropper(Y, X, ses, tau.sq = c("reml", "kNN"), H=1, 
       opt.method=c("MM", "optim"), control = list())
}
%- maybe also 'usage' for other objects documented here.
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
restricted maximum likelihood estimation (REML) or k-nearest neighbors.
}
  \item{H}{
The degree of the approximation used for the unbiased risk estimator. This must 
be an integer greater than or equal to 1. 
}
  \item{opt.method}{
The optimization method used for computing the regression coefficient estimates. 
This can be either an MM algorithm or the use of the \code{optim} function.       
}
  \item{control}{A list of control parameters specifying any changes to
      default values of the optimization procedures. Full names of control 
      list elements must be specified, otherwise, user-specifications are ignored.
	See *Details*.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x) 
{
  }
}
\keyword{ regression }
\keyword{ robust }
\concept{ empirical Bayes }
\concept{ loss-function driven estimation }
\concept{ unbiased risk estimation }
\concept{ ranking }
\concept{ mixed models }
