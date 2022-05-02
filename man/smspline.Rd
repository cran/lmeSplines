\name{smspline}
\alias{smspline}
\alias{smspline.v}
\title{Smoothing splines in NLME}
\description{Functions to generate matrices for a smoothing spline
  covariance structure, to enable fitting smoothing spline terms in LME/NLME. }
\usage{
  smspline(formula, data)
  smspline.v(time)
}
\arguments{
  \item{formula}{model formula with right hand side giving the spline covariate}
  \item{data}{optional data frame}
  \item{time}{spline `time' covariate to smooth over}
}
\details{A smoothing spline can be represented as a mixed model (Speed
  1991, Verbyla 1999).  The generated \eqn{Z}-matrix from \code{smspline()} 
  can be incorporated in the
  users's dataframe, then used in model formulae for LME random effects 
  terms at any level of grouping (see examples). The spline random terms
  are fitted in LME using an identity \code{'pdMat'} structure of the
  form \code{pdIdent(~Z - 1)}.
  The model formulation for a spline  in time (\code{t}) is as follows (Verbyla 1999):
\deqn{y = X_s \beta_s + Z_s u_s + e}
where
\eqn{X_s = [1 | t]} ,
\eqn{Z_s = Q(t(Q) Q)^{-1}} ,
and \eqn{u_s ~ N(0,G_s)}, is a set of random effects.
We transform the set of random effects \eqn{u_s} to independence with
\eqn{u_s = L v_s}, where \deqn{v_s ~ N(0,I \sigma^2_s)} is a set of
independent random effects. The \eqn{Z}-matrix is transformed
accordingly to \eqn{Z = Z_s L}, where \eqn{L} is the lower triangle of
the Choleski decomposition of \eqn{G_s}.

The function \code{smspline.v()} is called by \code{smspline()}, and can
be used to access the matrices \eqn{X_s, Q, G_s}. See Verbyla (1999)
for further information.

} 
\value{For \code{smspline()}, a \eqn{Z}-matrix with the same number of
rows as the
  data frame. After fitting, the LME model output gives a standard
  deviation parameter for the random effects, estimating
  \eqn{\sigma_s}. The smoothing parameter from the penalised likelihood
  formulation is \deqn{\lambda = \sigma^2/\sigma^2_s} % (check).
  
  For \code{smspline.v()}, a list of the form
  \item{Xs}{\eqn{X}-matrix for fixed effects part of the model}
  \item{Zs}{\eqn{Z}-matrix for random effects part of the model}
  \item{Q,Gs,R}{Matrices \eqn{Q, G_s, R} associated to the mixed-model form of
  the smoothing spline.}
}

\note{The time points for the smoothing spline basis are, by default,
  the unique values of the time covariate. This is the easiest approach,
  and model predictions at the fitted data points, can be obtained using
  \code{\link{predict.lme}}. By interpolation, using \code{\link{approx.Z}}, 
  the \eqn{Z}-matrix can be obtained for any set of time points and can 
  be used for fitting and/or prediction. (See examples).
Synopsis:data$Z <- smspline(formula1, data); fit <-lme(formula2, data, random= \dots)
}
\references{
  The correspondence between penalized likelihood formulations of
  smoothing splines and mixed models was pointed out by Speed (1991).
  The formulation used here for the mixed smoothing spline matrices are 
  given in Verbyla (1999). LME/NLME modelling is introduced in Pinheiro 
  and Bates (2000).

  Pinheiro, J. and Bates, D. (2000) Mixed-Effects Models in S and S-PLUS
  Springer-Verlag, New York.

  Speed, T. (1991) Discussion of ``That BLUP is a good thing: the
  estimation of random effects'' by G. Robinson.  Statist. Sci., 6, 42--44.

  Verbyla, A. (1999)  Mixed Models for Practitioners, Biometrics SA, Adelaide.

}

\author{Rod Ball \email{rod.ball@scionresearch.com} 
\url{https://www.scionresearch.com/}}
\seealso{\code{\link{approx.Z}} \code{\link{nlme}}}
\examples{
# smoothing spline curve fit
data(smSplineEx1)
# variable `all' for top level grouping
smSplineEx1$all <- rep(1,nrow(smSplineEx1))
# setup spline Z-matrix
smSplineEx1$Zt <- smspline(~ time, data=smSplineEx1)
fit1s <- lme(y ~ time, data=smSplineEx1,
    random=list(all=pdIdent(~Zt - 1)))
summary(fit1s)
plot(smSplineEx1$time,smSplineEx1$y,pch="o",type="n",
     main="Spline fits: lme(y ~ time, random=list(all=pdIdent(~Zt-1)))",
     xlab="time",ylab="y")
points(smSplineEx1$time,smSplineEx1$y,col=1)
lines(smSplineEx1$time, smSplineEx1$y.true,col=1)
lines(smSplineEx1$time, fitted(fit1s),col=2)

# fit model with cut down number of spline points
times20 <- seq(1,100,length=20)
Zt20 <- smspline(times20)
smSplineEx1$Zt20 <- approx.Z(Zt20,times20,smSplineEx1$time)
fit1s20 <- lme(y ~ time, data=smSplineEx1,
    random=list(all=pdIdent(~Zt20 - 1)))
# note: virtually identical df, loglik.
anova(fit1s,fit1s20)
summary(fit1s20)

# model predictions on a finer grid
times200 <- seq(1,100,by=0.5)
pred.df <- data.frame(all=rep(1,length(times200)),time=times200)
pred.df$Zt20 <- approx.Z(Zt20, times20,times200)
yp20.200 <- predict(fit1s20,newdata=pred.df)
lines(times200,yp20.200+0.02,col=4)


# mixed model spline terms at multiple levels of grouping
data(Spruce)
Spruce$Zday <- smspline(~ days, data=Spruce)
Spruce$all <- rep(1,nrow(Spruce))
# overall spline term, random plot and Tree effects
spruce.fit1 <- lme(logSize ~ days, data=Spruce,
                   random=list(all= pdIdent(~Zday -1),
                     plot=~1, Tree=~1))
# try overall spline term plus plot level linear + spline term
spruce.fit2 <- lme(logSize ~ days, data=Spruce,
                   random=list(all= pdIdent(~Zday - 1),
                     plot= pdBlocked(list(~ days,pdIdent(~Zday - 1))),
                     Tree = ~1))
anova(spruce.fit1,spruce.fit2)
summary(spruce.fit1)
}  
\keyword{smooth}
\keyword{models}
\keyword{regression}

