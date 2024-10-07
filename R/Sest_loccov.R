# Computes S-estimates by using the 
# fast S algorithm for multivariate location/covariance estimation
# implemented in the rrcov package 
# INPUT:
#   Y : response matrix (n x m)
#   bdp : breakdown point (<= 0.5)
# OUTPUT:
#   res$Mu : estimate of regression coefficients (or location vector)
#   res$Gamma : estimate of shape matrix 
#   res$scale : estimate of scale

#---------------------------------------------------------------------
Sest_loccov <- function(Y, bdp=.5, control=Scontrol(...), ...) {

    ## Computes Tukey's biweight psi function with constant c, divided by x, for all values in x
    scaledpsibiweight <- function(x, c){
        
        hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
        psi <- hulp*(abs(x)<c)
        
        return(psi)
    }

    ## performs vec-operation (stacks colums of a matrix into column-vector)
    vecop <- function(mat) {
    
        nr <- nrow(mat)
        nc <- ncol(mat)
        
        vecmat <- rep(0,nr*nc)
        for (col in 1:nc) {
            startindex <- (col-1)*nr+1
            vecmat[startindex:(startindex+nr-1)] <- mat[,col]
        }
        return(vecmat)
    }

    # reconstructs vecop'd matrix
    reconvec <- function(vec,ncol) {
    
        lcol <- length(vec)/ncol
        rec <- matrix(0,lcol,ncol)
        for (i in 1:ncol)
            rec[,i] <- vec[((i-1)*lcol+1):(i*lcol)]
        
        return(rec)
    }

    # converts control list into an rrcov control object for S-estimator
    convert <- function(control) {
        
        control = rrcov::CovControlSest(eps = control$convTol, maxiter = control$maxIt,
                nsamp = control$nsamp, trace = FALSE, method= "sfast")
        return(control)
    }

    Y <- as.matrix(Y)
    n <- nrow(Y)
    m <- ncol(Y)
    
    if(n <= m) 
        stop("number of observations too small (should have n > m)")
    
    Sests <- rrcov::CovSest(Y, bdp=bdp, control=convert(control))
    w <- scaledpsibiweight(sqrt(rrcov::getDistance(Sests)),Sests@cc)
    
    list(Mu=t(rrcov::getCenter(Sests)), 
         Gamma=rrcov::getShape(Sests), 
         scale=rrcov::getDet(Sests)^(1/(2*m)), 
         Sigma=rrcov::getCov(Sests), 
         c=Sests@cc, b=Sests@kp, w=w, outFlag=(!rrcov::getFlag(Sests)))
}


