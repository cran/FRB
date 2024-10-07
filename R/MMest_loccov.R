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
# --------------------------------------------------------------------

MMest_loccov <- function(Y, control=MMcontrol(...), ...) {

    ## Computes scaled Tukey's biweight psi function with constant c, divided by x, for all values in x
    scaledpsibiweight <- function(x, c) {
        hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
        psi <- hulp*(abs(x)<c)    
        return(psi)
    }
    
    ## converts control list into an rrcov control object for MM-estimator
    convert <- function(control) {
        control.S = rrcov::CovControlSest(eps = control$fastScontrols$convTol, 
        maxiter = control$fastScontrols$maxIt,nsamp = control$fastScontrols$nsamp, 
        trace = FALSE, method= "sfast")
        
        control = rrcov::CovControlMMest(bdp=control$bdp,eff=control$eff,sest=control.S,
        tolSolve = control$convTol.MM, maxiter = control$maxIt.MM)
        return(control)
    }

    Y <- as.matrix(Y)
    n <- nrow(Y)
    m <- ncol(Y)
    
    if(n <= m) 
        stop("number of observations too small (should have n > m)")
    
    MMests = rrcov::CovMMest(Y,eff.shape=control$shapeEff,control=convert(control))
    w <- scaledpsibiweight(sqrt(rrcov::getDistance(MMests)),MMests@c1)
    
    list(Mu=t(rrcov::getCenter(MMests)), 
        Gamma=rrcov::getShape(MMests), 
        scale=rrcov::getDet(MMests)^(1/(2*m)), 
        Sigma=rrcov::getCov(MMests), 
        c1=MMests@c1, 
        SMu=t(rrcov::getCenter(MMests@sest)),
        SGamma=rrcov::getShape(MMests@sest),
        SSigma=rrcov::getCov(MMests@sest),c0=MMests@sest@cc,
        b=MMests@sest@kp, w=w, outFlag=(!rrcov::getFlag(MMests)))
}


