
MMest_multireg <- function(X, Y, control=MMcontrol(...), ...)
{
# computes multivariate regression (M)M-estimates with auxiliary S-scale 
# INPUT:
# 	Y : response matrix (n x m)
# 	X : covariates matrix (n x p), possibly including intercept column
# 		(X = as.matrix(rep(1,n)) in case of location estimation)
# OUTPUT:
# 	res$Beta : MM-estimate of regression coefficients
#	  res$Gamma : MM-estimate of shape matrix 
#	  res$Sigma : MM-estimate of covariance matrix
# 	res$SBeta : S-estimate of regression coefficients
#   res$SGamma : S-estimate of shape matrix
#	  res$SSigma : S-estimate of covariance matrix
#   res$scale : S-estimate of scale
#
# breakdown point (bdp) set at .5
# efficiency : 95% efficiency for regression coefficients

# calls: Sest_multireg()
# --------------------------------------------------------------------
# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes scaled Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------
# (taken from Claudia Becker's Sstart0 program)

"chi.int" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (0,c1) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))

# --------------------------------------------------------------------

"loceff.bw" <- function(p, c1)
{  
# called by csolve.bw.MM(); computes efficiency corresponding to c1

alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))   
beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2

return( beta1^2 / alpha1 )

}

# --------------------------------------------------------------------

csolve.bw.MM <- function(p, eff)
{
# constant for second Tukey Biweight rho-function for MM, for fixed location-efficiency 

maxit <- 1000
eps <- 10^(-8)
diff <- 10^6
#ctest <- csolve.bw.asymp(p,.5)
ctest <- -.4024 + 2.2539 * sqrt(p) # very precise approximation for c corresponding to 50% bdp
iter <- 1
while ((diff>eps) & (iter<maxit)) {
    cold <- ctest
    ctest <- cold * eff / loceff.bw(p,cold)
    diff <- abs(cold-ctest)
    iter <- iter+1
}
return(ctest)

}

# --------------------------------------------------------------------
# -                       main function                              -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(Y)
m <- ncol(Y)
p <- ncol(X)

eff <- control$eff
bdp <- control$bdp
shapeEff <- control$shapeEff
fastScontrols <- control$fastScontrols
maxiter <- control$maxIt.MM
mtol <- control$convTol.MM

c1 <- csolve.bw.MM(m, eff)

# first compute 50% breakdown S-estimator
Sresult <- Sest_multireg(X, Y, bdp, fastScontrols)

auxscale <- Sresult$scale
newG <- Sresult$Gamma
newBeta <- Sresult$Beta
newR <- Y - X %*% newBeta
psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
newobj <- mean(rhobiweight(psres/auxscale,c1))
origobj <- newobj

# compute M-estimate with auxilliary scale through IRLS steps, starting
# from S-estimate
iteration <- 1
oldobj <- newobj + 1
while (((oldobj - newobj) > mtol) & (iteration < maxiter)) {
    oldobj <- newobj
    w <- scaledpsibiweight(psres/auxscale,c1)
    wbig <- matrix(rep(w,p),ncol=p)	
    wX <- X * wbig	
#    newBeta <- solve(t(wX) %*% X)%*%(t(wX) %*% Y)
    newBeta <- solve(crossprod(wX, X), crossprod(wX, Y))
    newG <- cov.wt(newR, wt=w, center=FALSE)$cov
    newG <- det(newG)^(-1/m)*newG
    newR <- Y - X %*% newBeta
    psres <- sqrt(mahalanobis(newR, rep(0,m), newG))
    newobj <- mean(rhobiweight(psres/auxscale,c1))
    iteration <- iteration+1
}

if (newobj <= origobj) {
    resBeta <- newBeta
    resshape <- newG
    rescovariance <- newG*auxscale^2
} else { # isn't supposed to happen
    resBeta <- Sresult$Beta
    resshape <- Sresult$Gamma
    rescovariance <- Sresult$Gamma*auxscale^2
}

return(list(Beta=resBeta, Gamma=resshape, Sigma=rescovariance, SBeta=Sresult$Beta,
        SGamma=Sresult$Gamma, SSigma=Sresult$Gamma*auxscale^2, scale=auxscale,
        c1 = c1, c0 = Sresult$c, b = Sresult$b))
 
}
