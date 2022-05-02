# lmeSplines
 lmeSplines R package (Version: 1.1-11) available on CRAN.

To install from GitHub:

`library("devtools"); install_github("agalecki/lmeSplines")`

Run example from `smspline` documentation:

````
library(lmeSplines)
data(smSplineEx1)

# variable `all` for top level grouping
smSplineEx1$all <- rep(1,nrow(smSplineEx1))
# setup spline Z-matrix
smSplineEx1$Zt <- smspline(~ time, data=smSplineEx1)
fit1s <- lme(y ~ time, data=smSplineEx1,
    random=list(all=pdIdent(~Zt - 1)))
summary(fit1s)

````
