plotDiag <- function(x, onepage = TRUE,...) {

currentAsk <- devAskNewPage(ask = NULL)
FRBres <- x

Y <- FRBres$Y
X <- FRBres$X
n <- nrow(X)
q <- ncol(Y)
# compute the residual distances:
resids <- Y - X %*% FRBres$Beta
if (!is.null(FRBres$int)) {   # for GS-estimates
    resids <- resids - as.matrix(rep(1,n))%*%FRBres$int 
    X <- cbind(rep(1,n), X)
}
residsD <- sqrt(mahalanobis(resids, rep(0,q), FRBres$Sigma))

# now perform multivariate location/scatter estimation via the same estimator
interceptdetection <- apply(X==1, 2, all)
withoutind <- (1:ncol(X))[interceptdetection==FALSE]
XwI <- X[,withoutind]

if (FRBres$method$est=="S")
  ests <- Sest_multireg(as.matrix(rep(1,n)), XwI, bdp=FRBres$method$bdp, control=FRBres$control)
else if (FRBres$method$est=="GS")
  ests <- GSest_multireg(as.matrix(rep(1,n)), XwI, bdp=FRBres$method$bdp, control=FRBres$control)
else
  ests <- MMest_multireg(as.matrix(rep(1,n)), XwI, control=FRBres$control)

Xloc <- ests$Beta
Xscatter <- ests$Sigma
leverageD <- sqrt(mahalanobis(XwI, Xloc, Xscatter))

# Least squares for comparison
B <- solve(crossprod(X), crossprod(X,Y))
R <- Y - X %*% B
S <- crossprod(R) /(n-(ncol(X)*q))
ClasresidsD <- sqrt(mahalanobis(R, rep(0,q), S))
ClasleverageD <- sqrt(mahalanobis(XwI, apply(XwI,2,mean), cov(XwI)))

##########################################

if (onepage) par(mfrow = c(2,1))  else   par(mfrow = c(1,1))

plot(leverageD, residsD, xlab="Robust distance in X-space", ylab="Robust residual distance", cex.lab=1.3, pch=20)
abline(h=sqrt(qchisq(.975, q)))
abline(v=sqrt(qchisq(.975, ncol(XwI))))

# identify outliers
inds <- 1:n
badinds <- inds[residsD > sqrt(qchisq(.975, q))]
if (length(badinds) <= 20) {
  for (j in badinds)
    text(leverageD[j], residsD[j], j, pos=4)
}

if (!onepage) devAskNewPage(ask = TRUE) 
      
plot(ClasleverageD, ClasresidsD, xlab="Mahalanobis distance in X-space", ylab="Residual distance", cex.lab=1.3, pch=20)
abline(h=sqrt(qchisq(.975, q)))
abline(v=sqrt(qchisq(.975, ncol(XwI))))
if (length(badinds) <= 20) {
  for (j in badinds)
    text(ClasleverageD[j], ClasresidsD[j], j, pos=4)
}

devAskNewPage(ask = currentAsk) 
}

