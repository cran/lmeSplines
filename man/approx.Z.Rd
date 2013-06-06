\title{Interpolating in smoothing spline Z-matrix columns}
\name{approx.Z}
\alias{approx.Z}
\description{
  Function to interpolate the \eqn{Z}-matrix for LME smoothing spline fits
  from one set of values of the \code{time} covariate to another.
}
\usage{
approx.Z(Z, oldtimes, newtimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{\eqn{Z}-matrix with rows corresponding to the sorted unique
    values of the \code{time} covariate }
  \item{oldtimes}{original (sorted) values for \code{time} covariate,
  corresponding to the rows of \eqn{Z}}
  \item{newtimes}{new (sorted) values for \code{time} covariate}
}
\details{
  Uses linear interpolation of each column of the \eqn{Z}-matrix, regarded as
  a function of time, with the times given by \code{oldtimes}.
}
\value{
  A matrix corresponding to an interpolated spline matrix. Can be used
  with \code{smspline()} for fitting LME splines with random effects
  corresponding to a different set of values of the \code{time} covariate 
  than those represented in the data, or as part of the   `newdata' argument 
  prediction from an existing model, to obtain predictions at points
  not represented in the data using \code{predict.lme()}.
}
\author{Rod Ball \email{rod.ball@scionresearch.com} 
  \url{www.scionresearch.com}}
\note{Linear interpolation works well here because the spline basis
  functions are approximately piecewise linear.}
\seealso{\code{\link{smspline}} \code{\link{nlme}}  }
\examples{
times1 <- 1:10
Zt1 <- smspline(~ times1)
times2 <- seq(1,10,by=0.1)
Zt2 <- approx.Z(Zt1,oldtimes=times1,newtimes=times2)
}
\keyword{manip}

