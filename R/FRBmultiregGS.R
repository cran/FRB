
FRBmultiregGS <- function(X,...) UseMethod("FRBmultiregGS")

FRBmultiregGS.formula <- function(formula, data, ...)
{
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    Y <- model.response(mf)

    z <- FRBmultiregGS.default(X, Y, ...)
    return(z)
}                                                 

FRBmultiregGS.default <- function(X, Y, R=999, bdp=0.5, conf=0.95, control=GScontrol(...), ...)
{
# performs multivariate regression based on multivariate GS estimates, with
# fast and robust bootstrap
#
# calls: GSest_multireg.R(), GSboot_multireg.R(), GSeinfs_multireg()
#
# INPUT :
# 	Y : response matrix (n x q)
# 	X : covariates matrix (n x p), possibly including intercept column
# 		(X = as.matrix(rep(1,n)) in case of location estimation)
#   R : number of bootstrap samples
#   bdp : breakdown point of GS-estimate (determines tuning parameters)
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$GSest : (list) result of GSmulti.R
#   res$GSboot : (list) result of multiGSregboot.R
#   res$Beta : (p*q) GS-Beta estimate
#   res$Sigma : (q*q) GS-Sigma estimate
#   res$SE : (p*q) bootstrap standard errors for each element of GS-Beta
#   res$CI.bca.lower : (p*q) lower bounds of 95% BCa intervals for each 
#                        element of GS-Beta
#   res$CI.bca.upper : (p*q) upper bounds of 95% BCa intervals for each 
#                        element of GS-Beta
#   res$CI.basic.lower : (p*q) lower bounds of 95% basic bootstrap intervals
#                          for each element of GS-Beta
#   res$CI.basic.upper : (p*q) upper bounds of 95% basic bootstrap intervals
#                          for each element of GS-Beta
                                                           

# --------------------------------------------------------------------

vecop <- function(mat) {
# performs vec-operation (stacks colums of a matrix into column-vector)

nr <- nrow(mat)
nc <- ncol(mat)

vecmat <- rep(0,nr*nc)
for (col in 1:nc) {
    startindex <- (col-1)*nr+1
    vecmat[startindex:(startindex+nr-1)] <- mat[,col]
}
return(vecmat)
}

# --------------------------------------------------------------------

reconvec <- function(vec,ncol) {
# reconstructs vecop'd matrix

lcol <- length(vec)/ncol
rec <- matrix(0,lcol,ncol)
for (i in 1:ncol)
    rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]

return(rec)
}

# --------------------------------------------------------------------
# -                        main function                             -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(Y)
q <- ncol(Y)
p <- ncol(X)

interceptdetection <- apply(X==1, 2, all)
#interceptind <- (1:p)[interceptdetection==TRUE]
#colnames(X)[interceptind] <- "(intercept)"
zonderint <- (1:p)[interceptdetection==FALSE]
Xzonderint <- X[,zonderint]
X <- Xzonderint
p<-ncol(X)

dimens <- p*q + q*q

if (is.null(colnames(Y)))
    colnames(Y) <- paste("Y",1:q,sep="")
if (is.null(colnames(X)))
    colnames(X) <- paste("X",1:p,sep="")

GSests <- GSest_multireg(X, Y, bdp=bdp, control=control)

GSBetawith <- GSests$Beta
GSint <- t(as.matrix(GSBetawith[1,]))
GSBeta <- as.matrix(GSBetawith[2:(p+1),])
GSSigma <- GSests$Sigma

if (q==1)
  {colnames(GSBeta) <- colnames(Y)
   colnames(GSBetawith) <- colnames(Y)
   colnames(GSSigma) <- colnames(Y)
   rownames(GSSigma) <- colnames(Y)
}

if (R<2) warning("argument R should be at least 2 to perform bootstrap inference; FRB is now skipped")

if (R>1) {
  bootres <- GSboot_multireg(X, Y, R=R,conf=conf, ests=GSests)

  stdsBeta <- reconvec(bootres$SE[1:(p*q)],q)

  #gives lower and upper bounds of BCA and basis bootstrap confidence intervals
  lowerlimitsBeta.bca <- reconvec(bootres$CI.bca[1:(p*q),1], q)
  upperlimitsBeta.bca <- reconvec(bootres$CI.bca[1:(p*q),2], q)
  lowerlimitsBeta.basic <- reconvec(bootres$CI.basic[1:(p*q),1], q)
  upperlimitsBeta.basic <- reconvec(bootres$CI.basic[1:(p*q),2], q)
}
else {
  bootres <- NULL
  stdsBeta <- NULL

  lowerlimitsBeta.bca <- NULL
  upperlimitsBeta.bca <- NULL
  lowerlimitsBeta.basic <- NULL
  upperlimitsBeta.basic <- NULL
}

####################################################################################

#method <- paste("Multivariate regression based on multivariate GS-estimates (breakdown point = ", bdp, ")", sep="")
method <- list(est="GS", bdp=bdp)

z <- list(est=GSests,bootest=bootres,Beta=GSBeta,intercept=GSint,Sigma=GSSigma, SE=stdsBeta, 
       CI.bca.lower=lowerlimitsBeta.bca,CI.bca.upper=upperlimitsBeta.bca,CI.basic.lower=lowerlimitsBeta.basic, CI.basic.upper=upperlimitsBeta.basic,
       conf=conf, method=method, control=control, X=X, Y=Y)

class(z) <- "FRBmultireg"
  
return(z)

}



