Sboot_multireg<- function(X, Y, R, conf=0.95, ests=Sest_multireg(X, Y))
{
# robust bootstrap for multivariate S-regression
# INPUT:
#   Y : n x q response matrix
#   X : n x p covariates matrix (ones(n,1) for location/shape estimation)
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
#   ests : result of Sest_multireg()
#
# OUTPUT: 
#   res$centered: ((p*q + q*q) x R) centered recomputations S-estimates:
#                             - first p*q rows: regression coeffients      
#                             - next q*q rows : covariance matrix
#                   (all in vec-form, columns stacked on top of each other) 
#   
#   res$vecest: ((p*q + q*q) x 1) original S estimates in vec-form
#   res$SE: ((p*q + q*q) x 1) bootstrap standard errors for elements in res$Sest
#   res$CI.bca: ((p*q + q*q) x 2) BCa confidence limits for elements in res$Sest
#   res$CI.basic: ((p*q + q*q) x 2) basic bootstrap confidence limits for elements in res$Sest
# --------------------------------------------------------------------

rhobiweight <- function(x,c)
{
# Computes Tukey's biweight rho function with constant c for all values in x

hulp <- x^2/2 - x^4/(2*c^2) + x^6/(6*c^4)
rho <- hulp*(abs(x)<c) + c^2/6*(abs(x)>=c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder1 <- function(x,c)
{
# Computes Tukey's biweight psi function with constant c for all values in x

hulp <- x - 2*x^3/(c^2) + x^5/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

rhobiweightder2 <- function(x,c)
{
# Computes derivative of Tukey's biweight psi function with constant c for all values in x

hulp <- 1 - 6*x^2/(c^2) + 5*x^4/(c^4)
rho <- hulp*(abs(x)<c)

return(rho)
}

# --------------------------------------------------------------------

commut <- function(p,m) {

# computes commutation matrix
# p = no. of rows
# m = no. of columns (of matrix which follows the commut matrix)

kompm <- matrix(0,p*m,p*m)
for (k in 1:(p*m)) {
    l <- (k - 1 - (ceiling(k/m)-1) * m) * p + ceiling(k/m)
    kompm[k,l] <- 1
}

return(kompm)

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

#------------------------------------------------------------------------------#

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
dimens <- p*q + q*q

Iq <- diag(rep(1,q))
Ip <- diag(rep(1,p))

c <- ests$c
b <- ests$b
Sigma0 <- ests$Sigma
Beta0 <- ests$Beta

Sinv <- solve(Sigma0)

###############################################################
###                 calculate jacobian                      ###
###############################################################

#           p*q          q*q        
#       -------------------------
#       |           |            |
#  p*q  | g1_Beta   | g1_Sigma   |
#       |           |            |
#       -------------------------
#       |           |            |
#  q*q  | g2_Beta   | g2_Sigma   |
#       |           |            |  
#       -------------------------

resmatrix <- Y - X %*% Beta0
divec <- sqrt(mahalanobis(resmatrix, rep(0,q), Sigma0))
divec[divec < 1e-5] <- 1e-5
udivec <- rhobiweightder1(divec,c)/divec
wdivec <- (rhobiweightder2(divec,c)*divec - rhobiweightder1(divec,c))/divec^3   
zdivec <- rhobiweightder2(divec,c)
wwdivec <- rhobiweightder1(divec,c)*divec - rhobiweight(divec,c)

#A <- t(X) %*% (matrix(rep(udivec,p),ncol=p) * X)
#B <- t(X) %*% (matrix(rep(udivec,q),ncol=q) * Y)
uX <- matrix(rep(udivec,p),ncol=p) * X
A <- crossprod(uX, X)
B <- crossprod(uX, Y)

term1a <- matrix(0, p*p,p*q)
term1b <- matrix(0,p*q,p*q)
term2a <- matrix(0,q*q,p*q)
term2b <- matrix(0,q*q,p*q)
term2c <- matrix(0,1,p*q)
term3a <- matrix(0,p*p,q*q)
term3b <- matrix(0,p*q,q*q)
term4a <- matrix(0,q*q,q*q)
term4b <- matrix(0,1,q*q)

for (i in 1:n) {
    Xi <- as.matrix(X[i,])
    Yi <- as.matrix(Y[i,])
    resi <- as.matrix(resmatrix[i,])
#    vecXiXi <- vecop(Xi %*% t(Xi))
#    vecXiYi <- vecop(Xi %*% t(Yi))
    vecXiXi <- vecop(tcrossprod(Xi))
    vecXiYi <- vecop(tcrossprod(Xi,Yi))
    vecresiresi <- vecop(tcrossprod(resi))
    wdi <- wdivec[i]
    zdi <- zdivec[i]
    udi <- udivec[i]
    veci <- vecop(Xi %*% t(resi) %*% Sinv)
    vecSi <- vecop(Sinv %*% resi %*% t(resi) %*% Sinv)

    term1a <- term1a + wdi * tcrossprod(vecXiXi, veci)
    term1b <- term1b + wdi * tcrossprod(vecXiYi, veci)

    term2a <- term2a + wdi * tcrossprod(vecresiresi, veci)
    term2b <- term2b + udi * ((kronecker(Iq,resi) + kronecker(resi,Iq)) %*% (kronecker(t(Xi),Iq) %*% commut(p,q))) 
    term2c <- term2c + zdi * t(veci)
    
    term3a <- term3a + wdi * tcrossprod(vecXiXi, vecSi)
    term3b <- term3b + wdi * tcrossprod(vecXiYi, vecSi)

    term4a <- term4a + wdi * tcrossprod(vecresiresi, vecSi)
    term4b <- term4b + zdi * t(vecSi)
}

Ainv <- solve(A)

partder1 <- (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% term1a) - (kronecker(Iq,Ainv) %*% term1b)   
partder2 <- -q/(b*n) * (term2a + term2b) + 1/(b*n) * vecop(Sigma0) %*% term2c
partder3 <- 1/2 * (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% term3a) - 1/2 * (kronecker(Iq,Ainv) %*% term3b)  
partder4 <- -q/(2*b*n) * term4a + 1/(2*b*n) * vecop(Sigma0) %*% term4b - 1/(b*n) * sum(wwdivec) * diag(rep(1,q*q))

jacobian <- cbind(rbind(partder1, partder2), rbind(partder3, partder4))

Idim <- diag(rep(1,dimens))
lincorrectmat <- solve(Idim-jacobian)

######################################################################

# put all estimates (coefs and covariances) in one column 
vecestim <- rep(0,dimens)
vecestim[1:(p*q)] <- vecop(Beta0)
vecestim[(p*q+1):dimens] <- vecop(Sigma0)

# to draw bootstrap samples 
#set.seed(2)
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens,R)  

for (r in 1:R) {
    Yst <- Y[bootmatrix[,r],]
    Xst <- X[bootmatrix[,r],]
    resmatrixst <- Yst - Xst %*% Beta0
    divecst <- sqrt(mahalanobis(resmatrixst, rep(0,q), Sigma0))
    divecst[divecst<1e-5] <- 1e-5
    udivecst <- rhobiweightder1(divecst,c)/divecst
    wwdivecst <- rhobiweightder1(divecst,c)*divecst - rhobiweight(divecst,c)     
    
    uXst <- matrix(rep(udivecst,p),ncol=p) * Xst
    Bst <- solve(crossprod(uXst, Xst), crossprod(uXst, Yst))
    uresst <- matrix(rep(udivecst,q),ncol=q) * resmatrixst
    Vst_term1 <- 1/(b*n) * q * crossprod(uresst, resmatrixst)
    Vst_term2 <- 1/(b*n)* sum(wwdivecst) * Sigma0
    Vst <- Vst_term1 - Vst_term2
    
    # list uncorrected bootstrap recomputations
    vecfst <- rep(0,dimens)
    vecfst[1:(p*q)] <- vecop(Bst)
    vecfst[(p*q+1):dimens] <- vecop(Vst)
    
    # compute centered, corrected fast bootstrap estimates
    fstbias <- vecfst - vecestim
    bootbiasmat[,r] <- lincorrectmat %*% fstbias  
}

#############################################################################

# compute bootstrap estimates of variance
Svariances <- apply(bootbiasmat, 1, var)

# sort bootstrap recalculations for constructing intervals
sortedSest <- t(apply(bootbiasmat, 1, sort))

# empirical inlfuences for computing a in BCa intervals, based on IF(S)
Einf <- Seinfs_multireg(X, Y, ests=ests)
inflE <- cbind(Einf$BetaS, Einf$covS)

normquan <- qnorm(1 - (1 - conf)/2)
estCIbca <- matrix(0,dimens,2)
estCIbasic <- matrix(0,dimens,2)
for (i in 1:(dimens)) {
  nofless <- length(sortedSest[i,sortedSest[i,]<=0])
  w <- qnorm(nofless/(R+1))
  a <- 1/6 * sum(inflE[,i]^3) / (sum(inflE[,i]^2)^(3/2))
  alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
  alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
  indexlow <- max((R+1)*alphatildelow,1)
  indexlow <- min(indexlow,R)
  indexhigh <- min((R+1)*alphatildehigh,R)
  indexhigh <- max(indexhigh,1)
  estCIbca[i,1] <- sortedSest[i,round(indexlow)] + vecestim[i]
  estCIbca[i,2] <- sortedSest[i,round(indexhigh)] + vecestim[i]
}

indexlow <- floor((1 - (1 - conf)/2) * R)
indexhigh <- ceiling((1 - conf)/2 * R)
estCIbasic[,1] <- vecestim - sortedSest[,indexlow]
estCIbasic[,2] <- vecestim - sortedSest[,indexhigh]

#############################################################################

return(list(centered=bootbiasmat, vecest=vecestim, SE=sqrt(Svariances), CI.bca=estCIbca, CI.basic=estCIbasic))

}

