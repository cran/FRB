
Sest_multireg <- function(X, Y, bdp=.5, control=Scontrol(...), ...)
{
# fast S algorithm for multivariate regression estimation
# INPUT:
# 	Y : response matrix (n x m)
# 	X : covariates matrix (n x p), possibly including intercept column
# 		(X = as.matrix(rep(1,n)) in case of location estimation)
# 	nsamp : number of elemental starts, e.g. 20 
# 	bdp : breakdown point (<= 0.5)
# OUTPUT:
# 	res$Beta : estimate of regression coefficients
#	  res$Gamma : estimate of shape matrix 
#	  res$scale : estimate of scale
#   res$Sigma : estimate of covariance matrix
# -------------------------------------------------------------------

IRLSstep <- function(X,Y,initialBeta, initialGamma, initialscale, k, c0, b, convTol)
{  
#convTol <- 1e-10
n <- nrow(X)
p <- ncol(X)
m <- ncol(Y)

Beta <- initialBeta
Res <- Y-X%*%Beta
psres <- sqrt(mahalanobis(Res, rep(0,m), initialGamma))
if (initialscale > 0)
    scale <- initialscale
else
    scale <- median(psres)/.6745

iter <- 0
betadiff <- 1

while ( (betadiff > convTol) & (iter < k) ) {
    iter <- iter + 1
    scale <- sqrt(scale^2 * mean(rhobiweight(psres/scale,c0))/b)
    w <- scaledpsibiweight(psres/scale,c0)
#    wbig <- diag(w)	
#    newBeta <- solve(t(X)%*%wbig%*%X)%*%(t(X)%*%wbig%*%Y)
    wbig <- matrix(rep(w,p),ncol=p)	
    wX <- X * wbig	
#    newBeta <- solve(t(wX) %*% X)%*%(t(wX) %*% Y)
    newBeta <- solve(crossprod(wX, X), crossprod(wX, Y))
    newGamma <- cov.wt(Res, wt=w, center=FALSE)$cov
    newGamma <- det(newGamma)^(-1/m)*newGamma
    Res <- Y-X%*%newBeta
    betadiff <- sum((newBeta-Beta)^2)/sum(Beta^2) # use 'sum' as a kind of norm
    Beta <- newBeta
    psres <- sqrt(mahalanobis(Res, rep(0,m), newGamma))
}
return(list( Beta = newBeta, Gamma = newGamma, scale = scale ))
}

#--------------------------------------------------------------------------  

scale1 <- function(u, b, c0, initialsc) 
{
# from Kristel's fastSreg
if (initialsc==0)
    initialsc = median(abs(u))/.6745
maxit <- 100
sc <- initialsc
i <- 0 
eps <- 1e-10
err <- 1
while  (( i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean(rhobiweight(u/sc,c0)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
}

return(sc)

}

#---------------------------------------------------------------------------------------

randomset <- function(tot,nel) {
ranset <- rep(0,nel)
for (j in 1:nel) {
 	num <- ceiling(runif(1)*tot)
   	if (j > 1) {
     		while (any(ranset==num)) 
     			num <- ceiling(runif(1)*tot)
	}
	ranset[j] <- num
}
return(ranset)
}

# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

psibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

# --------------------------------------------------------------------

scaledpsibiweight <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 2*x^2/(c^2) + x^4/(c^4)
psi <- hulp*(abs(x)<c)

return(psi)
}

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

#------------------------------------------------------------------------------#
#		     function needed to determine constant in biweight             #
#------------------------------------------------------------------------------#
# (taken from Claudia Becker's Sstart0 program)

"chi.int" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (0,c1) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * pchisq(c1^2, p + a))

"chi.int.p" <- function(p, a, c1)
  return(exp(lgamma((p + a)/2) - lgamma(p/2)) * 2^{a/2} * dchisq(c1^2, p + a) * 2 * c1)

"chi.int2" <- function(p, a, c1)
return(exp(lgamma((p + a)       
        #   partial expectation d in (c1,\infty) of d^a under chi-squared p
  /2) - lgamma(p/2)) * 2^{a/2} * (1 - pchisq(c1^2, p + a)))

"chi.int2.p" <- function(p, a, c1)
  return( - exp(lgamma((p + a)/2) - lgamma(p/2)) * 2^{a/2} * dchisq(c1^2, p + a) * 2 * c1)


"csolve.bw.asymp" <- function(p, r)
{
#   find biweight c that gives a specified breakdown
        c1 <- 9
        iter <- 1
        crit <- 100
        eps <- 0.00001
        while((crit > eps) & (iter < 100)) {
                c1.old <- c1
                fc <- erho.bw(p, c1) - (c1^2 * r)/6
                fcp <- erho.bw.p(p, c1) - (c1 * r)/3
                c1 <- c1 - fc/fcp
                if(c1 < 0)
                        c1 <- c1.old/2
                crit <- abs(fc) 
                iter <- iter + 1
        }
        return(c1)
}

"erho.bw" <- function(p, c1)
  return(chi.int(p,       #   expectation of rho(d) under chi-squared p
         2, c1)/2 - chi.int(p, 4, c1)/(2 * c1^2) + chi.int(p, 6, c1)/(6 * c1^4) + (
        c1^2 * chi.int2(p, 0, c1))/6)

"erho.bw.p" <- function(p, c1)
  return(chi.int.p(p,     #   derivative of erho.bw wrt c1
    2, c1)/2 - chi.int.p(p, 4, c1)/(2 * c1^2) + (2 * chi.int(p, 4, c1))/(2 * 
        c1^3) + chi.int.p(p, 6, c1)/(6 * c1^4) - (4 * chi.int(p, 6, c1))/(
        6 * c1^5) + (c1^2 * chi.int2.p(p, 0, c1))/6 + (2 * c1 * chi.int2(p,
        0, c1))/6)


#--------------------------------------------------------------------------#
#	       		end 'constant determining'				   #
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------------  
#--------------------------------------------------------------------------  
#--------------------------------------------------------------------------  

#set.seed(10)   # seed can be set, but be careful in simulations then...

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
m <- ncol(Y)

nsamp <- control$nsamp
bestr <- control$bestr # number of best solutions to keep for further C-steps
k <- control$k # number of C-steps on elemental starts
convTol <- control$convTol
maxIt <- control$maxIt

loop <- 1
c0 <- csolve.bw.asymp(m,bdp)
b <- erho.bw(m, c0)

bestbetas <- matrix(0, p*m, bestr)
bestgammas <- matrix(0, m*m, bestr)
bestscales <- 1e20 * rep(1,bestr)
sworst <- 1e20

while (loop <= nsamp) {
    # find a (p+m)-subset in general position.
    rankR <- 0
    itertest <- 0
    while ((rankR < m) && (itertest<200)) {
        ranset <- randomset(n,p+m)
        Xj <- X[ranset,]
        Yj <- Y[ranset,]
    	  Bj <- solve(crossprod(Xj), crossprod(Xj,Yj))
        Rj <- Yj - Xj %*% Bj
#	  qrXj <- qr(Xj)
	      qrRj <- qr(Rj)
        rankR <- qrRj$rank
	      itertest <- itertest + 1
    }
    if (itertest==200) stop("too many degenerate subsamples")

    Sj <- crossprod(Rj) /(p+m-1)
    Gj <- det(Sj)^(-1/m) * Sj
    # perform k steps of IRLS on elemental start
    res <- IRLSstep(X, Y, Bj, Gj, 0, k, c0, b, convTol)
  
    Betarw <- res$Beta
    Gammarw <- res$Gamma
    scalerw <- res$scale
    psresrw <- sqrt(mahalanobis(Y-X%*%Betarw, rep(0,m), Gammarw))
   if (loop > 1) {
        # check whether new Beta and new Gamma belong to the top best Betas; if so keep
        # Beta and Gamma with corresponding scale.
        if (mean(rhobiweight(psresrw/sworst,c0)) < b) {
            ss <- sort(bestscales, index.return=TRUE)
            ind <- ss$ix[bestr]
            bestscales[ind] <- scale1(psresrw, b, c0, scalerw)
            bestbetas[,ind] <- vecop(Betarw)
            bestgammas[,ind] <- vecop(Gammarw)
            sworst <- max(bestscales)
        }
    }
    else {
        bestscales[bestr] <- scale1(psresrw, b, c0, scalerw)
        bestbetas[,bestr] <- vecop(Betarw)
        bestgammas[,bestr] <- vecop(Gammarw)
    }
    loop <- loop + 1
}

ibest <- which.min(bestscales)
superbestscale <- bestscales[ibest]
superbestbeta <- reconvec(bestbetas[,ibest],m)
superbestgamma <- reconvec(bestgammas[,ibest],m)

# perform C-steps on best 'bestr' solutions, until convergence (or maximum maxIt steps) 
for (i in bestr:1) { 
    tmp <- IRLSstep(X, Y, reconvec(bestbetas[,i],m), reconvec(bestgammas[,i],m), bestscales[i], maxIt, c0, b, convTol)
    if (tmp$scale < superbestscale) {
        superbestscale <- tmp$scale;
        superbestbeta <- tmp$Beta;
        superbestgamma <- tmp$Gamma;
    }
}

return(list( Beta = superbestbeta, Gamma = superbestgamma, scale = superbestscale, 
                    Sigma = superbestgamma*superbestscale^2, b=b, c=c0))

}

