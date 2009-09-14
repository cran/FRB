
MMest_loccov <- function(Y, control=MMcontrol(...), ...)
{
# computes multivariate location and shape (M)M-estimates with auxiliary S-scale 
# INPUT:
#   Y = data
# OUTPUT:
#   res$Mu : MM-location estimate
#   res$Gamma : MM-shape matrix
#   res$Sigma : MM-covariance
#   res$SMu : S-location estimate
#   res$SGamma : S-shape matrix
#   res$SSigma : S-shape matrix
#   res$scale : S-scale
#
# calls: Sest_loccov()
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

"sigma1.bw" <- function(p, c1)
{  
# called by csolve.bw.MM()

gamma1.1 <- chi.int(p,2,c1) - 6*chi.int(p,4,c1)/(c1^2) + 5*chi.int(p,6,c1)/(c1^4)   
gamma1.2 <- chi.int(p,2,c1) - 2*chi.int(p,4,c1)/(c1^2) + chi.int(p,6,c1)/(c1^4)
gamma1 <- ( gamma1.1 + (p+1)*gamma1.2 ) / (p+2)

sigma1.0 <- chi.int(p,4,c1) - 4*chi.int(p,6,c1)/(c1^2) + 6*chi.int(p,8,c1)/(c1^4) - 4*chi.int(p,10,c1)/(c1^6) + chi.int(p,12,c1)/(c1^8)
return( sigma1.0 / (gamma1^2) * p/(p+2) )

}

# --------------------------------------------------------------------

"loceff.bw" <- function(p, c1)
{  
# called by csolve.bw.MM(); computes location efficiency corresponding to c1

alpha1 <- 1/p * (chi.int(p,2,c1) - 4*chi.int(p,4,c1)/(c1^2) + 6*chi.int(p,6,c1)/(c1^4) - 4*chi.int(p,8,c1)/(c1^6) + chi.int(p,10,c1)/(c1^8))   
beta1.1 <- chi.int(p,0,c1) - 2*chi.int(p,2,c1)/(c1^2) + chi.int(p,4,c1)/(c1^4)
beta1.2 <- chi.int(p,0,c1) - 6*chi.int(p,2,c1)/(c1^2) + 5*chi.int(p,4,c1)/(c1^4)
beta1 <- (1-1/p)*beta1.1 + 1/p*beta1.2

return( beta1^2 / alpha1 )

}

# --------------------------------------------------------------------

csolve.bw.MM <- function(p, eff, shape = TRUE)
{
# constant for second Tukey Biweight rho-function for MM, for fixed shape-efficiency 

maxit <- 1000
eps <- 10^(-8)
diff <- 10^6
#ctest <- csolve.bw.asymp(p,.5)
ctest <- -.4024 + 2.2539 * sqrt(p) # very precise approximation for c corresponding to 50% bdp
iter <- 1
while ((diff>eps) & (iter<maxit)) {
    cold <- ctest
    if (shape == TRUE)
        ctest <- cold * eff * sigma1.bw(p,cold)
    else
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
n <- nrow(Y)
q <- ncol(Y)

eff <- control$eff
bdp <- control$bdp
shapeEff <- control$shapeEff
fastScontrols <- control$fastScontrols
maxiter <- control$maxIt.MM
mtol <- control$convTol.MM

c1 <- csolve.bw.MM(q, eff, shape=shapeEff)

X <- as.matrix(rep(1,n))
 
# first compute the S-estimator
Sresult <- Sest_loccov(Y, bdp, fastScontrols)

auxscale <- Sresult$scale
newG <- Sresult$Gamma
newMu <- Sresult$Mu
newR <- Y - matrix(rep(newMu,n), n, byrow=TRUE)
psres <- sqrt(mahalanobis(newR, rep(0,q), newG))
newobj <- mean(rhobiweight(psres/auxscale,c1))
origobj <- newobj

# compute M-estimate with auxiliary scale through IRLS steps, starting
# from S-estimate
iteration <- 1
oldobj <- newobj + 1
while (((oldobj - newobj) > mtol) & (iteration < maxiter)) {
    oldobj <- newobj
    w <- scaledpsibiweight(psres/auxscale,c1)
    sqrtw <- sqrt(w)
    newMu <- crossprod(w, Y) / as.vector(crossprod(sqrtw))
    wbig <- matrix(rep(sqrtw,q),ncol=q) 
    wR <- newR * wbig   
    newG <- crossprod(wR)
    newG <- det(newG)^(-1/q)*newG
    newR <- Y - matrix(rep(newMu,n), n, byrow=TRUE)
    psres <- sqrt(mahalanobis(newR, rep(0,q), newG))
    newobj <- mean(rhobiweight(psres/auxscale,c1))
    iteration <- iteration+1
}

if (newobj <= origobj) {
    resloc <- newMu
    resshape <- newG
    rescovariance <- newG*auxscale^2
} else { # isn't supposed to happen
   resloc <- Sresult$Mu
   resshape <- Sresult$Gamma
   rescovariance <- Sresult$Gamma*auxscale^2
}

psres <- sqrt(mahalanobis(Y, resloc, rescovariance))
w <- scaledpsibiweight(psres,c1)
outFlag <- (psres > sqrt(qchisq(.975, q)))

return(list(Mu=resloc, Gamma=resshape, Sigma=rescovariance, SMu=Sresult$Mu,
        SGamma=Sresult$Gamma, SSigma=Sresult$Gamma*auxscale^2, scale=auxscale,
        c0=Sresult$c, b=Sresult$b, c1=c1, w=w, outFlag=outFlag))

}
