\name{smSplineEx1}
\alias{smSplineEx1}
\docType{data}
\title{Simulated data about a smooth curve}
\description{
  Simulated data to demonstrate smoothing spline curve fitting with
  \code{\link{smspline}} and \code{\link{lme}}
}
\usage{data(smSplineEx1)}
\format{
  A data frame with 100 observations on the following 4 variables.
  \describe{
    \item{time}{time covariate}
    \item{y}{simulated response values}
    \item{y.true}{true response values}
  }
}
\details{
100 data points were simulated about the curve  
\eqn{y = 10 - 6*\exp(-4t/100)}, with iid normal random errors
with standard deviation \eqn{1}.
}

\examples{
data(smSplineEx1)
}
\keyword{datasets}
