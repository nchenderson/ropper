\name{PPER}
\alias{PPER}
\title{
Population posterior expected ranks
}
\description{
Computes estimates of population posterior expected ranks given an estimate of 
the regression coefficients and an estimate of the random effects variance.
}
\usage{
PPER(beta.hat, Y, X, ses, tau.sq)
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A K x 1 vector of population posterior expected rank estimates.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{ropper}
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
\keyword{ regression }
\keyword{ robust }
\concept{ empirical Bayes }
\concept{ ranking }
% Use only one concept per line.
