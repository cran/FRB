
FRBmultiregMM <- function(X,...) UseMethod("FRBmultiregMM")

FRBmultiregMM.formula <- function(formula, data, ...)
{
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    Y <- model.response(mf)

    z <- FRBmultiregMM.default(X, Y, int = FALSE, ...)
    return(z)
}                                                 

FRBmultiregMM.default <- function(X, Y, int = TRUE, R=999, conf=0.95, control=MMcontrol(...), ...)
{
# performs multivariate regression based on multivariate MM estimates, with
# fast and robust bootstrap
#
# calls: MMest_multireg(), MMboot_multireg(), MMeinfs_multireg()
#
# INPUT :
# 	Y : response matrix (n x q)
# 	X : covariates matrix (n x p) or (n x (p-1))
#   int : logical; if TRUE, an intercept column is added
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
# OUTPUT :
#   res$est : (list) result of MMest_multireg()
#   res$bootest : (list) result of MMboot_multireg()
#   res$Beta : (p x q) MM estimate of the regression coefficient matrix
#   res$Sigma : (q x q) MM estimate of the error covariance matrix
#   res$SE : (p x q) bootstrap standard errors for MM-estimate Beta
#   res$CI.bca.lower : (p x q) 95% BCa lower limits for Beta
#   res$CI.bca.upper : (p x q) 95% BCa upper limits for Beta
#   res$CI.basic.lower : (p x q) 95% basic bootstrap lower limits for Beta
#   res$CI.basic.upper : (p x q) 95% basic bootstrap upper limits for Beta

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
#bdp <- .5

if (is.null(colnames(Y)))
    colnames(Y) <- paste("Y",1:q,sep="")
if (is.null(colnames(X)))
    colnames(X) <- paste("X",1:p,sep="")

interceptdetection <- apply(X==1, 2, all)
interceptind <- (1:p)[interceptdetection==TRUE]
colnames(X)[interceptind] <- "(intercept)"

if (!any(interceptdetection) & int){
    X <- cbind(rep(1,n),X)
    p <- p + 1    
    colnames(X)[1] <- "(intercept)"
}
dimens <- p*q + q*q

MMests <- MMest_multireg(X, Y, control=control)
MMBeta <- MMests$Beta
MMSigma <- MMests$Sigma

if (R<2) warning("argument R should be at least 2 to perform bootstrap inference; FRB is now skipped")

if (R>1) {
  bootres <- MMboot_multireg(X, Y, R=R, conf=conf, ests=MMests)

  stdsBeta <- reconvec(bootres$SE[1:(p*q)],q)

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

###################################################################################

#method <- paste("Multivariate regression based on multivariate MM-estimates (bdp = ", control$bdp, ", eff = ", control$eff, ")", sep="")

method <- list(est="MM", bdp=control$bdp, eff=control$eff)

z <- list(est=MMests, bootest=bootres, Beta=MMBeta, Sigma=MMSigma, SE=stdsBeta, CI.bca.lower=lowerlimitsBeta.bca, 
        CI.bca.upper=upperlimitsBeta.bca, CI.basic.lower=lowerlimitsBeta.basic, CI.basic.upper=upperlimitsBeta.basic,
        conf=conf, method=method, control=control, X=X, Y=Y)
  
class(z) <- "FRBmultireg"
  
return(z)      

}