# Copyright (C) 2003 Roderick D. Ball
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# Or see www.gnu.org/copyleft/gpl.html

#' Smoothing Splines in NLME
#'
#' @description
#' Functions to generate matrices for a smoothing spline covariance structure,
#' enabling the fitting of smoothing spline terms in linear mixed-effects models
#' (LME) or nonlinear mixed-effects models (NLME). A smoothing spline can be
#' represented as a mixed model, as described by Speed (1991) and Verbyla (1999).
#' The generated Z-matrix can be incorporated into a data frame and used in LME
#' random effects terms with an identity covariance structure
#' (\code{pdIdent(~Z - 1)}).
#'
#' The model formulation for a spline in time (\code{t}) is:
#' \deqn{y = X_s \beta_s + Z_s u_s + e}
#' where \eqn{X_s = [1 | t]}, \eqn{Z_s = Q (t(Q) Q)^{-1}}, and
#' \eqn{u_s \sim N(0, G_s)} is a set of random effects. The random effects are
#' transformed to independence via \eqn{u_s = L v_s}, where
#' \eqn{v_s \sim N(0, I \sigma^2_s)} and \eqn{L} is the lower triangle of the
#' Cholesky decomposition of \eqn{G_s}. The Z-matrix is transformed to
#' \eqn{Z = Z_s L}.
#'
#' @param formula Model formula with the right-hand side specifying the spline
#'   covariate (e.g., \code{~ time}). Must contain exactly one variable.
#' @param data Optional data frame containing the variable specified in
#'   \code{formula}. If not provided, the formula is evaluated in the current
#'   environment.
#' @param time Numeric vector of spline time covariate values to smooth over.
#'
#' @return
#' For \code{smspline}, a Z-matrix with the same number of rows as the input data
#' frame or vector, representing the random effects design matrix for the
#' smoothing spline. After fitting an LME model, the standard deviation parameter
#' for the random effects estimates \eqn{\sigma_s}, and the smoothing parameter
#' is \eqn{\lambda = \sigma^2 / \sigma^2_s}.
#'
#' For \code{smspline.v}, a list containing:
#' \describe{
#'   \item{Xs}{Matrix for fixed effects, with columns \code{[1 | t]}.}
#'   \item{Zs}{Matrix for random effects, computed as \code{Q (t(Q) \%*\% Q)^-1 L}.}
#'   \item{Q}{Matrix used in the spline formulation.}
#'   \item{Gs}{Covariance matrix for the random effects.}
#'   \item{R}{Cholesky factor of \code{Gs}.}
#' }
#'
#' @author Rod Ball <rod.ball@scionresearch.com>
#' @references
#' Pinheiro, J. and Bates, D. (2000) \emph{Mixed-Effects Models in S and S-PLUS}.
#' Springer-Verlag, New York.
#'
#' Speed, T. (1991) Discussion of "That BLUP is a good thing: the estimation of
#' random effects" by G. Robinson. \emph{Statistical Science}, 6, 42--44.
#'
#' Verbyla, A. (1999) \emph{Mixed Models for Practitioners}. Biometrics SA,
#' Adelaide.
#'
#' @note
#' The time points for the smoothing spline basis are, by default, the unique
#' values of the time covariate. Model predictions at the fitted data points can
#' be obtained using \code{predict.lme}. For predictions at different time points,
#' use \code{\link{approx.Z}} to interpolate the Z-matrix.
#'
#' @seealso \code{\link{approx.Z}}, \code{\link[nlme]{lme}}
#' @examples
#' # Smoothing spline curve fit
#' data(smSplineEx1)
#' smSplineEx1$all <- rep(1, nrow(smSplineEx1))
#' smSplineEx1$Zt <- smspline(~ time, data = smSplineEx1)
#' fit1s <- lme(y ~ time, data = smSplineEx1,
#'              random = list(all = pdIdent(~ Zt - 1)))
#' summary(fit1s)
#' plot(smSplineEx1$time, smSplineEx1$y, pch = "o", type = "n",
#'      main = "Spline fits: lme(y ~ time, random = list(all = pdIdent(~ Zt - 1)))",
#'      xlab = "time", ylab = "y")
#' points(smSplineEx1$time, smSplineEx1$y, col = 1)
#' lines(smSplineEx1$time, smSplineEx1$y.true, col = 1)
#' lines(smSplineEx1$time, fitted(fit1s), col = 2)
#'
#' # Fit model with reduced number of spline points
#' times20 <- seq(1, 100, length = 20)
#' Zt20 <- smspline(times20)
#' smSplineEx1$Zt20 <- approx.Z(Zt20, times20, smSplineEx1$time)
#' fit1s20 <- lme(y ~ time, data = smSplineEx1,
#'                random = list(all = pdIdent(~ Zt20 - 1)))
#' anova(fit1s, fit1s20)
#' summary(fit1s20)
#'
#' # Model predictions on a finer grid
#' times200 <- seq(1, 100, by = 0.5)
#' pred.df <- data.frame(all = rep(1, length(times200)), time = times200)
#' pred.df$Zt20 <- approx.Z(Zt20, times20, times200)
#' yp20.200 <- predict(fit1s20, newdata = pred.df)
#' lines(times200, yp20.200 + 0.02, col = 4)
#'
#' # Mixed model spline terms at multiple levels of grouping
#' data(Spruce)
#' Spruce$Zday <- smspline(~ days, data = Spruce)
#' Spruce$all <- rep(1, nrow(Spruce))
#' spruce.fit1 <- lme(logSize ~ days, data = Spruce,
#'                    random = list(all = pdIdent(~ Zday - 1),
#'                                  plot = ~ 1, Tree = ~ 1))
#' spruce.fit2 <- lme(logSize ~ days, data = Spruce,
#'                    random = list(all = pdIdent(~ Zday - 1),
#'                                  plot = pdBlocked(list(~ days, pdIdent(~ Zday - 1))),
#'                                  Tree = ~ 1))
#' anova(spruce.fit1, spruce.fit2)
#' summary(spruce.fit1)
#'
#' @importFrom nlme lme pdIdent pdBlocked
#' @importFrom stats model.frame
#' @export
#' @keywords smooth models regression
smspline <- function(formula, data) {
  # Copyright (C) 2003 Roderick D. Ball
  # generate the Z matrix
  if (is.vector(formula)) {
    x <- formula
  } else {
    if (missing(data)) {
      mf <- model.frame(formula)
    } else {
      mf <- model.frame(formula, data)
    }
    if (ncol(mf) != 1) stop("formula can have only one variable")
    x <- mf[, 1]
  }

  x.u <- sort(unique(x))
  Zx <- smspline.v(x.u)$Zs
  Zx[match(x, x.u), ]
}

#' @rdname smspline
#' @export
smspline.v <- function(time) {
  # Copyright (C) 2003 Roderick D. Ball
  ### generate the X and Z matrices for a mixed model
  ### formulation of a smoothing spline as per Verbyla's notes 5.3.2
  ### Verbyla, A.P. Mixed models for practitioners. University of Adelaide and
  ### The South Australian Research and Development Institute, 1999, 115pp.
  ### limitation: knot points identical to data
  ### smoothing penalty \lambda_s \int g''(t) dt
  ### lambda_s = sigma^2/sigma^2_s
  ### y = X_s beta_s + Z_s u_s + e
  ### X_s = [1 | t]
  ### Z_s = Q (t(Q) %*%Q)^-1
  ### u_s ~ N(0, G_s)
  ### y[i] = g(t[i]) + e[i]
  ### let u_s = L v_s
  ### so cov(u) = G_s = L L'
  t1 <- sort(unique(time))
  p <- length(t1)
  h <- diff(t1)
  h1 <- h[1:(p - 2)]
  h2 <- h[2:(p - 1)]
  Q <- matrix(0, nrow = p, ncol = p - 2)
  Q[cbind(1:(p - 2), 1:(p - 2))] <- 1/h1
  Q[cbind(1 + 1:(p - 2), 1:(p - 2))] <- -1/h1 - 1/h2
  Q[cbind(2 + 1:(p - 2), 1:(p - 2))] <- 1/h2
  Gs <- matrix(0, nrow = p - 2, ncol = p - 2)
  Gs[cbind(1:(p - 2), 1:(p - 2))] <- 1/3 * (h1 + h2)
  Gs[cbind(1 + 1:(p - 3), 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
  Gs[cbind(1:(p - 3), 1 + 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
  # Zus = Q %*% (t(Q)%*%Q)^-1
  Zus <- t(solve(t(Q) %*% Q, t(Q)))
  #R _ choleski(Gs)
  R <- chol(Gs, pivot = FALSE)
  tol <- max(1e-12, 1e-8 * mean(diag(R)))
  if (sum(abs(diag(R))) < tol)
    stop("singular G matrix")
  Zvs <- Zus %*% t(R)
  list(Xs = cbind(rep(1, p), t1), Zs = Zvs, Q = Q, Gs = Gs, R = R)
}

#' Interpolating in Smoothing Spline Z-Matrix Columns
#'
#' @description
#' Interpolates the Z-matrix for LME smoothing spline fits from one set of time
#' covariate values to another using linear interpolation of each column of the
#' Z-matrix, regarded as a function of time.
#'
#' @param Z Z-matrix with rows corresponding to the sorted unique values of the
#'   time covariate (e.g., from \code{\link{smspline}} or \code{\link{smspline.v}}).
#' @param oldtimes Numeric vector of original (sorted) time covariate values
#'   corresponding to the rows of \code{Z}.
#' @param newtimes Numeric vector of new time covariate values to interpolate to.
#'
#' @return A matrix with the same number of columns as \code{Z} and rows
#'   corresponding to \code{newtimes}, containing the interpolated Z-matrix values.
#'   This can be used with \code{\link{smspline}} for fitting LME splines with
#'   random effects at different time points or as part of the \code{newdata}
#'   argument in \code{\link[nlme]{predict.lme}} for predictions at new points.
#'
#' @author Rod Ball <rod.ball@scionresearch.com>
#' @note
#' Linear interpolation works well because the spline basis functions are
#' approximately piecewise linear.
#'
#' @seealso \code{\link{smspline}}, \code{\link[nlme]{lme}},
#' \code{\link[nlme]{predict.lme}}
#' @examples
#' times1 <- 1:10
#' Zt1 <- smspline(~ times1)
#' times2 <- seq(1, 10, by = 0.1)
#' Zt2 <- approx.Z(Zt1, oldtimes = times1, newtimes = times2)
#'
#' @importFrom stats approx
#' @export
#' @keywords manip
approx.Z <- function(Z, oldtimes, newtimes) {
  # Copyright (C) 2003 Roderick D. Ball
  # linear interpolation of Z matrix by column.
  oldt.u <- sort(unique(oldtimes))
  if (length(oldt.u) != length(oldtimes) ||
      any(oldt.u != oldtimes)) {
    Z <- Z[match(oldt.u, oldtimes), ]
    oldtimes <- oldt.u
  }
  apply(Z, 2, function(u, oldt, newt) {
    approx(oldt, u, xout = newt)$y
  }, oldt = oldtimes, newt = newtimes)
}

#' Simulated Data for Smoothing Spline Curve Fitting
#'
#' @description
#' Simulated dataset to demonstrate smoothing spline curve fitting with
#' \code{\link{smspline}} and \code{\link[nlme]{lme}}. The data consists of 100
#' observations simulated around the curve \eqn{y = 10 - 6 \exp(-4t/100)} with
#' independent normal random errors (standard deviation = 1).
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{time}{Time covariate.}
#'   \item{y}{Simulated response values.}
#'   \item{y.true}{True response values.}
#' }
#'
#' @examples
#' data(smSplineEx1)
#' str(smSplineEx1)
#'
#' @docType data
#' @keywords datasets
"smSplineEx1"
