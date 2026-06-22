\name{mlb2025}
\alias{mlb2025}
\docType{data}
\title{
Baseball Winning Percentages in 2025
}
\description{
Data on the number of team wins during the 2025 
major league baseball (MLB) season. The data includes 
the following additional team characteristics: payroll, strength of schedule,
and ballpark difficulty.
}
\usage{data("mlb2025")}
\format{
  A data frame with 29 observations on the following 8 variables.
  \describe{
    \item{\code{Team}}{a character vector containing team names.}
    \item{\code{Wins}}{a numeric vector giving the number of total wins over 
    the 2025 MLB season.}
    \item{\code{Losses}}{a numeric vector giving the number of total losses over
    the 2025 MLB season.}
    \item{\code{WinProp}}{a numeric vector with the proportion of wins for 
    each team.}
    \item{\code{StdErrSq}}{a numeric vector with the standard errors squared 
    associated with the win proportions.}
    \item{\code{Payroll}}{a numeric vector containing the projected payroll 
    (in millions) at the start of the 2025 baseball season.}
    \item{\code{SOS}}{a numeric vector measuring the 2025 strength of schedule.}
    \item{\code{ParkFactors}}{a numeric vector measuring the overall difficulty 
    of each team's home ballpark. Park difficulty is averaged over the 2024-2026
    seasons.}
  }
}
\details{
All MLB teams except the Athletics were included in this data.
}
\source{
Total wins, total losses, and 2025 strength of schedule were obtained 
from \url{https://baseball-reference.com}. The park factor index
was retrieved from \url{https://baseballsavant.mlb.com/}. 
Payroll projections for 2025 were obtained from \url{https://blogs.fangraphs.com/the-2025-payrolls-and-beyond/}.
}
\examples{
data(mlb2025)
str(mlb2025)

## Fit model with Payroll, SOS, ParkFactors as covariates
X1 <- model.matrix(~ Payroll + SOS + ParkFactors, 
                   data=mlb2025)
rank_mod <- ropper(Y = mlb2025$WinProp, 
                   X=X1, 
                   ses=mlb2025$StdErrSq)

## Look at estimated regression coefficients
rank_mod$coef

## Look at ROPPER percentiles estimates
## Rank teams according to best percentiles
mlb_estperc <- data.frame(Team=mlb2025$Team, 
                      percentile=round(rank_mod$pepp, 4))
mlb_estperc <- mlb_estperc[order(-mlb_estperc$percentile),]

head(mlb_estperc)
}
\keyword{datasets}
