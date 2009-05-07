MMboot_multireg <- function(X, Y, R, conf=0.95, ests = MMest_multireg(X, Y))
{
# robust bootstrap for multivariate MM regression 
# INPUT:
#   Y : n x q response matrix
#   X : n x p covariates matrix (ones(n,1) for location/shape estimation)
#   R : number of bootstrap samples
#   conf : confidence level for bootstrap intervals
#   ests : result of multiMM_regression
#
# OUTPUT: 
#   res$centered: (2*(p*q + q*q) x R) centered recomputations of MM and S-estimates:
#                             - first p*q rows: MM location (or regression)      
#                             - next q*q rows : MM shape matrix 
#                             - next q*q rows : S covariance matrix
#                             - final p*q rows: S location (or regression)
#                   (all in vec-form, columns stacked on top of each other) 
#   res$vecest: (2*(p*q + q*q) x 1) original estimates in vec-form
#   res$SE: (2*(p*q + q*q) x 1) bootstrap standard errors for elements in res$vecest
#   res$CI.bca: (2*(p*q + q*q) x 2) BCa confidence limits for elements in res$vecest
#   res$CI.basic: (2*(p*q + q*q) x 2) basic bootstrap confidence limits for elements in res$vecest

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
# Computes Tukey's biweight psi function with constant c for all values in x

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

# --------------------------------------------------------------------
# -                         main function                            -
# --------------------------------------------------------------------

Y <- as.matrix(Y)
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)
dimens <- p*q + q*q

Iq <- diag(rep(1,q))
Ip <- diag(rep(1,p))

c0 <- ests$c0
b <- ests$b
c1 <- ests$c1

MMBeta <- ests$Beta
MMSigma <- ests$Sigma
Beta0 <- ests$SBeta
Sigma0 <- ests$SSigma

MMSinv <- solve(MMSigma)
S0inv <- solve(Sigma0)
auxscalesq <- det(Sigma0)^(1/q)
auxscale <- sqrt(auxscalesq)
MMGamma <- auxscalesq^(-1)*MMSigma
MMGinv <- solve(MMGamma)

##########################################################################
###                        calculate jacobian                          ###
##########################################################################

#           p*q          q*q          q*q            p*q  
#       -----------------------------------------------------
#       |           |            |            |             |
#  p*q  | g1_MMBeta | g1_MMGamma | g1_Sigma_0 |      0      |
#       |           |            |            |             |
#       -----------------------------------------------------
#       |           |            |            |             |
#  q*q  | g2_MMBeta | g2_MMGamma | g2_Sigma_0 |      0      |
#       |           |            |            |             |  
#       -----------------------------------------------------
#       |           |            |            |             |
#  q*q  |     0     |      0     | g3_Sigma_0 |  g3_Beta_0  |
#       |           |            |            |             |
#       -----------------------------------------------------
#       |           |            |            |             |
#  p*q  |     0     |      0     | g4_Sigma_0 |  g4_Beta_0  |
#       |           |            |            |             |
#       -----------------------------------------------------
#             1            2            3            4 

# first fill up lower right part: the S-part : g3, g4

restildematrix <- Y - X %*% Beta0
ditildevec <- sqrt(mahalanobis(restildematrix, rep(0,q), Sigma0))
ditildevec[ditildevec < 1e-5] <- 1e-5
uditildevec <- rhobiweightder1(ditildevec,c0)/ditildevec
wditildevec <- (rhobiweightder2(ditildevec,c0)*ditildevec - rhobiweightder1(ditildevec,c0))/ditildevec^3   
zditildevec <- rhobiweightder2(ditildevec,c0)
wwditildevec <- rhobiweightder1(ditildevec,c0)*ditildevec - rhobiweight(ditildevec,c0)

utildeX <- matrix(rep(uditildevec,p),ncol=p) * X
Atilde <- crossprod(utildeX, X)
Btilde <- crossprod(utildeX, Y)

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
    resi <- as.matrix(restildematrix[i,])
    vecXiXi <- vecop(tcrossprod(Xi))
    vecXiYi <- vecop(tcrossprod(Xi,Yi))
    vecresiresi <- vecop(tcrossprod(resi))
    wdi <- wditildevec[i]
    zdi <- zditildevec[i]
    udi <- uditildevec[i]
    veci <- vecop(Xi %*% t(resi) %*% S0inv)
    vecSi <- vecop(S0inv %*% resi %*% t(resi) %*% S0inv)

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

Atildeinv <- solve(Atilde)

partder1 <- (t(kronecker(Btilde,Ip)) %*% kronecker(t(Atildeinv),Atildeinv) %*% term1a) - (kronecker(Iq,Atildeinv) %*% term1b)   
partder2 <- -q/(b*n) * (term2a + term2b) + 1/(b*n) * vecop(Sigma0) %*% term2c
partder3 <- 1/2 * (t(kronecker(Btilde,Ip)) %*% kronecker(t(Atildeinv),Atildeinv) %*% term3a) - 1/2 * (kronecker(Iq,Atildeinv) %*% term3b)  
partder4 <- -q/(2*b*n) * term4a + 1/(2*b*n) * vecop(Sigma0) %*% term4b - 1/(b*n) * sum(wwditildevec) * diag(rep(1,q*q))

Part33 <- partder4
Part34 <- partder2
Part43 <- partder3
Part44 <- partder1

# end S-part

# now g1, g2 

resmatrix <- Y - X %*% MMBeta
divec <- sqrt(mahalanobis(resmatrix, rep(0,q), MMGamma))
divec[divec < 1e-5] <- 1e-5
udivec <- rhobiweightder1(divec/auxscale,c1)/divec
wdivec <- (rhobiweightder2(divec/auxscale,c1)*divec/auxscale - rhobiweightder1(divec/auxscale,c1))/divec^3 
vdivec <- rhobiweightder2(divec/auxscale,c1)

uX <- matrix(rep(udivec,p),ncol=p) * X
A <- crossprod(uX, X)
B <- crossprod(uX, Y)
V <- t(resmatrix) %*% (matrix(rep(udivec,q),ncol=q) * resmatrix)

termg1MMBetaa <- matrix(0,p*p,p*q);
termg1MMBetab <- matrix(0,p*q,p*q);
termg1MMGammaa <- matrix(0,p*p,q*q);
termg1MMGammab <- matrix(0,p*q,q*q);
termg2MMBetaa <- matrix(0,q*q,p*q);
termg2MMBetab <- matrix(0,q*q,p*q);
termg2MMGamma <- matrix(0,q*q,q*q);
termg1MMSigma0a <- matrix(0,p*p,1);
termg1MMSigma0b <- matrix(0,p*q,1);
termg2MMSigma0 <- matrix(0,q*q,1);

for (i in 1:n) {
    Xi <- as.matrix(X[i,])
    Yi <- as.matrix(Y[i,])
    resi <- as.matrix(resmatrix[i,])
    vecXiXi <- vecop(tcrossprod(Xi))
    vecXiYi <- vecop(tcrossprod(Xi,Yi))
    vecresiresi <- vecop(tcrossprod(resi))
    wdi <- wdivec[i]
    udi <- udivec[i]
    vdi <- vdivec[i]
    veci <- vecop(Xi %*% t(resi) %*% MMGinv)
    vecSi <- vecop(MMGinv %*% resi %*% t(resi) %*% MMGinv)

    termg1MMBetaa <- termg1MMBetaa + wdi * tcrossprod(vecXiXi, veci)
    termg1MMBetab <- termg1MMBetab + wdi * tcrossprod(vecXiYi, veci)

    termg1MMGammaa <- termg1MMGammaa + wdi * tcrossprod(vecXiXi, vecSi)
    termg1MMGammab <- termg1MMGammab + wdi * tcrossprod(vecXiYi, vecSi)

    termg2MMBetaa <- termg2MMBetaa + wdi * tcrossprod(vecresiresi, veci)
    termg2MMBetab <- termg2MMBetab + udi * ((kronecker(Iq,resi) + kronecker(resi,Iq)) %*% (kronecker(t(Xi),Iq) %*% commut(p,q))) 

    termg2MMGamma <- termg2MMGamma + wdi * tcrossprod(vecresiresi, vecSi)
    
    termg1MMSigma0a <- termg1MMSigma0a + vdi * vecXiXi
    termg1MMSigma0b <- termg1MMSigma0b + vdi * vecXiYi

    termg2MMSigma0 <- termg2MMSigma0 + vdi * vecresiresi
}

Ainv <- solve(A)

Part11 <- (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% termg1MMBetaa) - (kronecker(Iq,Ainv) %*% termg1MMBetab)
Part12 <- 1/2*(t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% termg1MMGammaa) - 1/2 * (kronecker(Iq,Ainv) %*% termg1MMGammab)
Part21 <- -det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% (termg2MMBetaa + termg2MMBetab)
Part22 <- -1/2 * det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% termg2MMGamma
Part13 <- -1/2/q/auxscale * (t(kronecker(B,Ip)) %*% kronecker(t(Ainv),Ainv) %*% termg1MMSigma0a - kronecker(Iq,Ainv) %*% termg1MMSigma0b) %*% t(vecop(t(S0inv)))   
Part23 <- -1/2/q/auxscale * det(V)^(-1/q) * (diag(rep(1,q*q)) - 1/q * vecop(V) %*% t(vecop(t(solve(V))))) %*% termg2MMSigma0 %*% t(vecop(t(S0inv)))

Part14 <- matrix(0,p*q,p*q)
Part24 <- matrix(0,q*q,p*q)
Part31 <- matrix(0,q*q,p*q)
Part32 <- matrix(0,q*q,q*q)
Part41 <- matrix(0,p*q,p*q)
Part42 <- matrix(0,p*q,q*q)

col1 <- rbind(Part11, Part21, Part31, Part41)
col2 <- rbind(Part12, Part22, Part32, Part42)
col3 <- rbind(Part13, Part23, Part33, Part43)
col4 <- rbind(Part14, Part24, Part34, Part44)
jacobian <- cbind(col1, col2, col3, col4)

Idim <- diag(rep(1,dimens*2))
lincorrectmat <- solve(Idim-jacobian)

######################################################################

# put all estimates (coefs and covariances) in one column 
vecestim <- rep(0,dimens*2)
vecestim[1:(p*q)] <- vecop(MMBeta)
vecestim[(p*q+1):dimens] <- vecop(MMGamma)
vecestim[(dimens+1):(dimens+(q*q))] <- vecop(Sigma0)
vecestim[(dimens+q*q+1):(dimens*2)] <- vecop(Beta0)

# to draw bootstrap samples 
# set.seed(2)
bootmatrix <- matrix(sample(n,R*n,replace=TRUE),ncol=R)

bootbiasmat <- matrix(0,dimens*2,R)  

for (r in 1:R) {
    Yst <- Y[bootmatrix[,r],]
    Xst <- X[bootmatrix[,r],]
    resmatrixst <- Yst - Xst %*% MMBeta
    divecst <- sqrt(mahalanobis(resmatrixst, rep(0,q), MMSigma))
    divecst[divecst<1e-5] <- 1e-5
    udivecst <- rhobiweightder1(divecst,c1)/divecst
    
    restildematrixst <- Yst - Xst %*% Beta0 
    ditildevecst<- sqrt(mahalanobis(restildematrixst, rep(0,q), Sigma0))
    ditildevecst[ditildevecst<1e-5] <- 1e-5
    uditildevecst <- rhobiweightder1(ditildevecst,c0)/ditildevecst
    wwditildevecst <- rhobiweightder1(ditildevecst,c0)*ditildevecst - rhobiweight(ditildevecst,c0)     
    
    uXst <- matrix(rep(udivecst,p),ncol=p) * Xst
    Bst <- solve(crossprod(uXst, Xst), crossprod(uXst, Yst))
    uresst <- matrix(rep(udivecst,q),ncol=q) * resmatrixst
    Gst <- crossprod(uresst, resmatrixst)
    Gst <- det(Gst)^(-1/q) * Gst    
    Vst <- auxscalesq * Gst
    
    utilderesst <- matrix(rep(uditildevecst,q),ncol=q) * restildematrixst
    V0st_term1 <- 1/(b*n) * q * crossprod(utilderesst, restildematrixst)
    V0st_term2 <- 1/(b*n)* sum(wwditildevecst) * Sigma0
    V0st <- V0st_term1 - V0st_term2
    
    utildeXst <- matrix(rep(uditildevecst,p),ncol=p) * Xst
    B0st <- solve(crossprod(utildeXst, Xst), crossprod(utildeXst, Yst))
    
    # list uncorrected bootstrap recomputations
    vecfst <- rep(0,dimens*2)
    vecfst[1:(p*q)] <- vecop(Bst)
    vecfst[(p*q+1):dimens] <- vecop(Gst)
    vecfst[(dimens+1):(dimens+q*q)] <- vecop(V0st)
    vecfst[(dimens+q*q+1):(dimens*2)] <- vecop(B0st)
    
    # compute centered, corrected fast bootstrap estimates
    fstbias <- vecfst - vecestim
    bootbiasmat[,r] <- lincorrectmat %*% fstbias  
}

# compute bootstrap estimates of variance
MMvariances <- apply(bootbiasmat, 1, var)

# sort bootstrap recalculations for constructing intervals
sortedMMest <- t(apply(bootbiasmat, 1, sort))

# empirical inlfuences for computing a in BCa intervals, based on IF(MM)
Einf <- MMeinfs_multireg(X, Y, ests=ests)
inflE <- cbind(Einf$Beta, Einf$shape, Einf$covS, Einf$BetaS)  

normquan <- qnorm(1 - (1 - conf)/2)
estCIbca <- matrix(0,dimens*2,2)
estCIbasic <- matrix(0,dimens*2,2)
for (i in 1:(dimens)) {
  nofless <- length(sortedMMest[i,sortedMMest[i,]<=0])
  w <- qnorm(nofless/(R+1))
  a <- 1/6 * sum(inflE[,i]^3) / (sum(inflE[,i]^2)^(3/2))
  alphatildelow <- pnorm(w+(w-normquan)/(1-a*(w-normquan)))
  alphatildehigh <- pnorm(w+(w+normquan)/(1-a*(w+normquan)))
  indexlow <- max((R+1)*alphatildelow,1)
  indexlow <- min(indexlow,R)
  indexhigh <- min((R+1)*alphatildehigh,R)
  indexhigh <- max(indexhigh,1)
  estCIbca[i,1] <- sortedMMest[i,round(indexlow)] + vecestim[i]
  estCIbca[i,2] <- sortedMMest[i,round(indexhigh)] + vecestim[i]
}

indexlow <- floor((1 - (1 - conf)/2) * R)
indexhigh <- ceiling((1 - conf)/2 * R)
estCIbasic[,1] <- vecestim - sortedMMest[,indexlow]
estCIbasic[,2] <- vecestim - sortedMMest[,indexhigh]

#############################################################################

return(list(centered=bootbiasmat, vecest=vecestim, SE=sqrt(MMvariances), CI.bca=estCIbca, CI.basic=estCIbasic))

}


