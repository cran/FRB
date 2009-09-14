
Sest_loccov <- function(Y, bdp=.5, control=Scontrol(...),...)
{
# fast S algorithm for multivariate location/covariance estimation
# INPUT:
#   Y : response matrix (n x m)
#   bdp : breakdown point (<= 0.5)
# OUTPUT:
#   res$Mu : estimate of regression coefficients (or location vector)
#   res$Gamma : estimate of shape matrix 
#   res$scale : estimate of scale

#---------------------------------------------------------------------

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
#                function needed to determine constant in biweight             #
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

# -------------------------------------------------------------------

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
#                      end 'constant determining'                          #
#--------------------------------------------------------------------------#

IRLSstep <- function(Y, initialMu, initialGamma, initialscale, k, c0, b, convTol)
{  
#convTol <- 1e-10
n <- nrow(Y)
m <- ncol(Y)

Mu <- initialMu
Res <- Y - matrix(rep(Mu,n), n, byrow=TRUE)
psres <- sqrt(mahalanobis(Res, rep(0,m), initialGamma))
if (initialscale > 0)
    scale <- initialscale
else
    scale <- median(psres)/.6745

iter <- 0
mudiff <- 1

while ( (mudiff > convTol) & (iter < k) ) {
    iter <- iter + 1
    scale <- sqrt(scale^2 * mean(rhobiweight(psres/scale,c0))/b)
    w <- scaledpsibiweight(psres/scale,c0)
    sqrtw <- sqrt(w)
    if(qr(Res[w>0,])$rank < m) stop("Too many points on a hyperplane!")
    newMu <- crossprod(w, Y) / as.vector(crossprod(sqrtw))
    wbig <- matrix(rep(sqrtw,m),ncol=m) 
    wRes <- Res * wbig  
    newGamma <- crossprod(wRes)
    newGamma <- det(newGamma)^(-1/m)*newGamma
    Res <- Y - matrix(rep(newMu,n), n, byrow=TRUE)
    mudiff <- sum((newMu-Mu)^2)/sum(Mu^2) # use 'sum' as a kind of norm
    Mu <- newMu
    psres <- sqrt(mahalanobis(Res, rep(0,m), newGamma))
}
return(list( Mu = newMu, Gamma = newGamma, scale = scale ))
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

# --------------------------------------------------------------------
# -                       main function                              -
# --------------------------------------------------------------------

#set.seed(10) # seed can be set, but be careful in simulations then...

Y <- as.matrix(Y)
n <- nrow(Y)
m <- ncol(Y)

nsamp <- control$nsamp
bestr <- control$bestr # number of best solutions to keep for further C-steps
k <- control$k # number of C-steps on elemental starts
convTol <- control$convTol
maxIt <- control$maxIt

# standardize the data
medY <- apply(Y,2,median)
madY <- apply(Y,2,mad)
Y <- Y - matrix(rep(medY,n), n, byrow=TRUE)
Y <- Y / matrix(rep(madY,n), n, byrow=TRUE)

if (n<=m) stop("number of observations too small (should have n > m)")

loop <- 1
c0 <- csolve.bw.asymp(m,bdp)
b <- erho.bw(m, c0)

bestmus <- matrix(0, m, bestr)
bestgammas <- matrix(0, m*m, bestr)
bestscales <- 1e20 * rep(1,bestr)
sworst <- 1e20   # magic number!!

while (loop <= nsamp) {
    # find a (p+m)-subset in general position.
    # find a (m+1)-subset in general position.
    rankcov <- m-1
    itertest <- 0
    while ((rankcov < m) && (itertest<200)) {
        ranset <- randomset(n,m+1)
        Yj <- Y[ranset,]
        Sj <- cov(Yj)
        rankcov <- qr(Sj)$rank
        itertest <- itertest + 1
    }
    if (itertest==200) stop("too many degenerate subsamples")

    Muj <- apply(Yj, 2, mean)
    Gj <- det(Sj)^(-1/m)*Sj
    # perform k steps of IRLS on elemental start
    res <- IRLSstep(Y, Muj, Gj, 0, k, c0, b, convTol)
  
    Murw <- res$Mu
    Gammarw <- res$Gamma
    scalerw <- res$scale
    psresrw <- sqrt(mahalanobis(Y - matrix(rep(Murw,n), n, byrow=TRUE), rep(0,m), Gammarw))
    if (loop > 1) {
        # check whether new Mu and new Gamma belong to the top best Mus; if so keep
        # Mu and Gamma with corresponding scale.
        if (mean(rhobiweight(psresrw/sworst,c0)) < b) {
            ss <- sort(bestscales, index.return=TRUE)
            ind <- ss$ix[bestr]
            bestscales[ind] <- scale1(psresrw, b, c0, scalerw)
            bestmus[,ind] <- vecop(as.matrix(Murw))
            bestgammas[,ind] <- vecop(Gammarw)
            sworst <- max(bestscales)
        }
    }
    else {
        bestscales[bestr] <- scale1(psresrw, b, c0, scalerw)
        bestmus[,bestr] <- vecop(as.matrix(Murw))
        bestgammas[,bestr] <- vecop(Gammarw)
    }
    loop <- loop + 1
}

ibest <- which.min(bestscales)
superbestscale <- bestscales[ibest]
superbestmu <- reconvec(bestmus[,ibest],m)
superbestgamma <- reconvec(bestgammas[,ibest],m)

# perform C-steps on best 'bestr' solutions, until convergence (or maximum maxIt steps) 
for (i in bestr:1) { 
    tmp <- IRLSstep(Y, reconvec(bestmus[,i],m), reconvec(bestgammas[,i],m), bestscales[i], maxIt, c0, b, convTol)
    if (tmp$scale < superbestscale) {
        superbestscale <- tmp$scale;
        superbestmu <- tmp$Mu;
        superbestgamma <- tmp$Gamma;
    }
}

psres <- sqrt(mahalanobis(Y, superbestmu, superbestgamma))/superbestscale
w <- scaledpsibiweight(psres,c0)
outFlag <- (psres > sqrt(qchisq(.975, m)))

# Unstandardize the results
UstMu <- superbestmu * madY + medY
superbestsigma <- superbestgamma * superbestscale^2
UstSigma <- diag(madY) %*% superbestsigma %*% diag(madY)
dimnames(UstSigma) <- dimnames(superbestsigma) 
UstGamma <- det(UstSigma)^(-1/m)*UstSigma
Ustscale <- det(UstSigma)^(1/2/m)

return(list( Mu = UstMu, Gamma = UstGamma, scale = Ustscale, Sigma = UstGamma*Ustscale^2, c=c0, b=b, w=w, outFlag=outFlag))

}
