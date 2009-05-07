GSest_multireg <- function(X,Y, bdp=.5, control=GScontrol(...), ...)
{
# A fast procedure to compute a GS-estimator based on the algorithm
# proposed by Salibian-Barrera, M. and Yohai, V.J. (2005),
# "A fast algorithm for S-regression estimates".
# 
#  Input:
#       y : response matrix (n x q) 
#       x : covariates matrix (n x p), possibly including intercept column
#       nsamp : number of sub-samples (default=200)
#       bdp : breakdown point (default=0.5)
#       
# Default arguments (can be redefined at the beginning of the program)
#       k = number of refining iterations in ea.  subsample (def=2)
#       bestr = number of "best betas" to remember
#               from the subsamples. These will be later
#               iterated until convergence (def=5)
#
#   Output is a list with components
#       beta : slope estimate of regression coefficients (p,q)  
#       covariance : robust estimate of scatter (p,p)
#       scale : estimate of scale (1,1)
#---------------------------------------------------------------------------


resmultiGS<-function(x,y,initialbeta,initialgamma,k,b,cc,initialscale,convTol)  {  
    # does "k" IRWLS refining steps from "initial.beta"
    #
    # if "initial.scale" is present, it's used, o/w the MAD is used
    # k = number of refining steps
    # b and cc = tuning constants of the equation
    # 
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)

    beta <- initialbeta
    res <- y-x%*%beta
    places <- t(combn(1:n,2))
    ngroot <- nrow(places)
    term1resvector <- as.matrix(res[places[,1],])
    term2resvector <- as.matrix(res[places[,2],])
    diffres <- term1resvector-term2resvector
    rdis <- sqrt(mahalanobis(diffres, rep(0,q), initialgamma))
    term1xvector <- x[places[,1],]
    term2xvector <- x[places[,2],]
    diffx <- term1xvector-term2xvector
    term1yvector <- y[places[,1],]
    term2yvector <- y[places[,2],]
    diffy <- term1yvector-term2yvector
    if (initialscale > 0)
    scale <- initialscale
    else
    scale <- median(rdis)/.6745
    
    gamma <- initialgamma

    for (i in 1:k) {
        # do one step of the iterations to solve for the scale
        scale <- sqrt( scale^2 * mean( rhobiweight( rdis / scale, cc ) ) / b)
        
        # now do one step of IRWLS with the "improved scale"
        w <- scaledpsibiweight( rdis/scale, cc )
        
        if (p>0) {
          wbig <- matrix(rep(w,p),ncol=p)	
          wdiffx <- diffx * wbig
          newbeta <- solve(crossprod(wdiffx, diffx), crossprod(wdiffx, diffy))
        }
        else
          newbeta <- beta
          
        newgamma <- cov.wt(diffres, wt=w, center=FALSE)$cov
        newgamma <- det(newgamma)^(-1/q)*newgamma
 
        # check for convergence
       if (p>0) {  
          if (sum((newbeta-beta)^2)/sum(beta^2) < convTol)       break
       }
       else {
          if (sum(sum((newgamma-gamma)^2))/sum(sum(gamma^2)) < convTol)       break
       }
        res <- y-x%*%newbeta 
        term1resvector <- as.matrix(res[places[,1],])
        term2resvector <- as.matrix(res[places[,2],])
        diffres <- term1resvector-term2resvector
        rdis <- sqrt(mahalanobis(diffres, rep(0,q), newgamma))
        gamma <- newgamma
        beta <- newbeta
    }    
    
return(list( betarw = beta, gammarw = newgamma, scalerw = scale,rdis=rdis ))
}
    
#------------------------------------------------------------------------
scalemultiGS <- function(u, b, cc, initialsc)
{
# from Kristel's fastSreg

if (initialsc==0)
 { initialsc <- median(abs(u))/.6745}
 
# find the scale, full iterations
maxit <- 200
# magic number alert
sc <- initialsc
i <- 0
eps <- 1e-20
# magic number alert
err <- 1
while ((i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean( rhobiweight( u / sc, cc ) ) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
}

return(sc)
}
#--------------------------------------------------------------------------  


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

#-----------------------------------------------------------------------

erf <- function(x)
{
  uitk <- 2*pnorm(x*sqrt(2))-1 
  return(uitk)
}

#-----------------------------------------------------------------------
 TbscGS <- function(alpha,p)
{
# constant for Tukey Biweight S 

talpha <- sqrt(qchisq(1-alpha,p))
maxit <- 1000 
eps <- 1e-8
diff <- 1e6
ctest <- talpha
iter <- 1
while ((diff>eps) * (iter<maxit))
 {   cold <- ctest
    if (alpha >= 0.50)
        {ctest <- TbsbGS(cold,p)/(1-alpha^2)}
    #bdp kleiner dan 0.50
    else
        {ctest <- TbsbGS(cold,p)/(1-(1-alpha)^2)}
    diff <- abs(cold-ctest)
    iter <- iter+1
  }
 return(ctest)
}

#---------------------------------------------------------------------------------------------
TbsbGS <- function(c,p)
{

if (p==1)
{  y1 <- -1/3*(60*exp(-1/4*c^2)*c+exp(-1/4*c^2)*c^5-60*pi^(1/2)*erf(1/2*c)-3*pi^(1/2)*erf(1/2*c)*c^4-8*exp(-1/4*c^2)*c^3+18*c^2*pi^(1/2)*erf(1/2*c))/(c^4*pi^(1/2))
   y2 <- -1/6*c^2*erf(1/2*c)+1/6*c^2
   res <- (6/c)*(y1+y2)
}
else 
{   if (p==2)
      {
       tus <- -1/6*(-384+96*c^2-12*c^4+384*exp(-1/4*c^2)+exp(-1/4*c^2)*c^6)/c^4+1/6*c^2*exp(-1/4*c^2)
       res <- (6/c)*tus
      }

   else 
     {  if (p==3)
        {tus <- -1/6*(2*exp(-1/4*c^2)*c^5*pi-40*exp(-1/4*c^2)*c^3*pi+840*exp(-1/4*c^2)*c*pi-840*erf(1/2*c)*pi^(3/2)+180*erf(1/2*c)*pi^(3/2)*c^2+exp(-1/4*c^2)*c^7*pi-18*erf(1/2*c)*pi^(3/2)*c^4)/(pi^(3/2)*c^4)+1/6*c^2*(exp(-1/4*c^2)*c*pi^2-pi^(5/2)*erf(1/2*c)+pi^(5/2))/pi^(5/2)
         res <- (6/c)*tus
        }
        else 
        {    if (p==4)
             {tus <- -1/24*(-96*c^4-6144+1152*c^2+6144*exp(-1/4*c^2)+384*c^2*exp(-1/4*c^2)+4*exp(-1/4*c^2)*c^6+exp(-1/4*c^2)*c^8)/c^4+1/24*c^2*exp(-1/4*c^2)*(4+c^2);
             res <- (6/c)*tus
             }
             else 
             {    if (p==5)
                  {tus <- -1/36*(12*exp(-1/4*c^2)*c^5*pi^2+exp(-1/4*c^2)*c^9*pi^2+2520*erf(1/2*c)*pi^(5/2)*c^2+6*exp(-1/4*c^2)*c^7*pi^2-180*erf(1/2*c)*pi^(5/2)*c^4-15120*pi^(5/2)*erf(1/2*c)+15120*exp(-1/4*c^2)*c*pi^2)/(pi^(5/2)*c^4)+1/36*c^2*(pi^4*exp(-1/4*c^2)*c^3+6*pi^4*exp(-1/4*c^2)*c-6*pi^(9/2)*erf(1/2*c)+6*pi^(9/2))/pi^(9/2)
                  res <- (6/c)*tus
                  }
                  else 
                  {    if (p==6)
                      {tus <- -1/192*(-1152*c^4-122880+18432*c^2+384*exp(-1/4*c^2)*c^4+122880*exp(-1/4*c^2)+12288*exp(-1/4*c^2)*c^2+32*exp(-1/4*c^2)*c^6+8*exp(-1/4*c^2)*c^8+exp(-1/4*c^2)*c^10)/c^4+1/192*c^2*exp(-1/4*c^2)*(32+8*c^2+c^4)
                      res <- (6/c)*tus
                      }
                      else 
                     {    if (p==7)
                          {tus <- -1/360*(60*exp(-1/4*c^2)*c^7*pi^3+45360*erf(1/2*c)*pi^(7/2)*c^2+504*exp(-1/4*c^2)*c^5*pi^3+10*exp(-1/4*c^2)*c^9*pi^3+332640*exp(-1/4*c^2)*c*pi^3+10080*exp(-1/4*c^2)*c^3*pi^3-2520*erf(1/2*c)*pi^(7/2)*c^4+exp(-1/4*c^2)*c^11*pi^3-332640*erf(1/2*c)*pi^(7/2))/(pi^(7/2)*c^4)+1/360*c^2*(pi^6*exp(-1/4*c^2)*c^5+10*pi^6*exp(-1/4*c^2)*c^3+60*pi^6*exp(-1/4*c^2)*c-60*pi^(13/2)*erf(1/2*c)+60*pi^(13/2))/pi^(13/2)
                           res <- (6/c)*tus
                          }
                          else 
                         {     if (p==8)
                               {tus <- -1/2304*(-18432*c^4-2949120+368640*c^2+18432*exp(-1/4*c^2)*c^4+2949120*exp(-1/4*c^2)+368640*exp(-1/4*c^2)*c^2+768*exp(-1/4*c^2)*c^6+96*exp(-1/4*c^2)*c^8+exp(-1/4*c^2)*c^12+12*exp(-1/4*c^2)*c^10)/c^4+1/2304*c^2*exp(-1/4*c^2)*(384+96*c^2+12*c^4+c^6)
                                res <- (6/c)*tus
                               }
                               else 
                              {    if (p==9)
                                  {tus <- -1/5040*(23184*exp(-1/4*c^2)*c^5*pi^4+exp(-1/4*c^2)*c^13*pi^4+140*exp(-1/4*c^2)*c^9*pi^4+14*exp(-1/4*c^2)*c^11*pi^4+997920*erf(1/2*c)*pi^(9/2)*c^2-8648640*erf(1/2*c)*pi^(9/2)+1224*exp(-1/4*c^2)*c^7*pi^4+443520*exp(-1/4*c^2)*c^3*pi^4-45360*erf(1/2*c)*pi^(9/2)*c^4+8648640*exp(-1/4*c^2)*c*pi^4)/(pi^(9/2)*c^4)+1/5040*c^2*(pi^8*exp(-1/4*c^2)*c^7+14*pi^8*exp(-1/4*c^2)*c^5+140*pi^8*exp(-1/4*c^2)*c^3+840*pi^8*exp(-1/4*c^2)*c-840*pi^(17/2)*erf(1/2*c)+840*pi^(17/2))/pi^(17/2)
                                   res <-(6/c)*tus
                                  }
                                  else 
                                 {     if (p==10)
                                       {tus <- -1/36864*(-368640*c^4+8847360*c^2-82575360+737280*c^4*exp(-1/4*c^2)+11796480*c^2*exp(-1/4*c^2)+82575360*exp(-1/4*c^2)+30720*exp(-1/4*c^2)*c^6+1920*exp(-1/4*c^2)*c^8+16*exp(-1/4*c^2)*c^12+192*exp(-1/4*c^2)*c^10+exp(-1/4*c^2)*c^14)/c^4+1/36864*c^2*exp(-1/4*c^2)*(6144+1536*c^2+192*c^4+16*c^6+c^8)
                                        res <- (6/c)*tus
                                       }
else 
   {    if (p>10)
   {res <- TbsbGSbenader(c,p)

                                           
                                            }
                                            
                                            }
                                       }
                                  }
                             }
                        }
                    }
               }
         }
     }
}
}  


#------------------------------------------------------------------------------------------------------------------


TbsbGSbenader <- function(c,p) {

integralpart1 <- function(r) {((r^2)/2-(r^4)/(2*c^2)+(r^6)/(6*c^4))*exp((-r^2)/4)*r^(p-1)}

integralpart2<- function(r) {exp((-r^2)/4)*r^(p-1)}

tus1<-(2/(gamma(p/2)*2^p))*integrate(integralpart1,0,c)$value
tus2<-((c^2)/(gamma(p/2)*3*2^p))*integrate(integralpart2,c,Inf)$value
tustot<-tus1+tus2
res<-(6/c)*tustot
return(res)
}

#---------------------------------------------------------------------------
IRLSlocation <- function(xmat,covmat,bdp,cc)
{
n <- nrow(xmat)
p <- ncol(xmat)

neem <- sample(n,p+1)
xsub <- xmat[neem,]

initmu <- apply(xsub,2,mean)
initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))

initobj <- mean(rhobiweight(initrdis,cc))

weights <- scaledpsibiweight(initrdis,cc)

itertest <- 0
while ((sum(weights)==0) && (itertest<500)) {
    
    neem <- sample(n,p+1)
    xsub <- xmat[neem,]
   
    initmu <- apply(xsub,2,mean)
    initrdis <- sqrt(mahalanobis(xmat, initmu, covmat))
    initobj <- mean(rhobiweight(initrdis,cc))
    weights <- scaledpsibiweight(initrdis,cc)
    itertest <- itertest + 1
}
if (itertest==500) stop("could not find suitable starting point for IRLS for intercept")   
   
sqrtweights <- sqrt(weights)
munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))

rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))

objnieuw <- mean(rhobiweight(rdisnieuw,cc))

iter <- 0
while (((abs(initobj/objnieuw)-1) > 10^(-15)) && (iter < 100)) {
    initobj <- objnieuw
    initmu <- munieuw
    initrdis <- rdisnieuw
    weights <- scaledpsibiweight(initrdis,cc)
    sqrtweights <- sqrt(weights)
    munieuw <- crossprod(weights, xmat) / as.vector(crossprod(sqrtweights))
    rdisnieuw <-sqrt(mahalanobis(xmat,munieuw,covmat))
    objnieuw <- mean(rhobiweight(rdisnieuw,cc))
    iter <- iter + 1
}
return(munieuw)
}



# ---------------------------------------------------------------------------------------------------
# -                        main function                                                           -
# --------------------------------------------------------------------------------------------------

p<-ncol(X)
interceptdetection <- apply(X==1, 2, all)
zonderint <- (1:p)[interceptdetection==FALSE]
xzonderint <- X[,zonderint]
X <- xzonderint

nsamp <- control$nsamp
bestr <- control$bestr # number of best solutions to keep for further C-steps
k <- control$k # number of C-steps on elemental starts
convTol <- control$convTol
maxIt <- control$maxIt

X <- as.matrix(X)
Y <- as.matrix(Y)
n <- nrow(X)
ngroot <- choose(n,2)
p <- ncol(X)
q <- ncol(Y)
cc <- TbscGS(bdp,q)
b <- (cc/6)*TbsbGS(cc,q)

#set.seed(10)

bestbetas <- matrix(0,p*q,bestr)
bestgammas <- matrix(0,q*q,bestr)
bestscales <- 1e20 * rep(1,bestr)
sworst <- 1e20
xextra <- cbind(rep(1,n),X)

for(i in 1:nsamp)
  { # 
    # find a (p+q+1)-subset in general position.
    #
    
    rankR <- 0
    itertest <- 0

    while ((rankR < q) && (itertest<200)) {
	  ranset <- sample(n,p+q+1)
        xj <- xextra[ranset,]
        yj <- Y[ranset,]
    	  beta <- solve(crossprod(xj), crossprod(xj,yj))
        res <- yj - xj %*% beta
	      qrRj <- qr(res)
        rankR <- qrRj$rank
        itertest <- itertest + 1
    }
    if (itertest==200) stop("too many degenerate subsamples")
	
    Smat <- crossprod(res)/(p+q-1)
    Cmat <- det(Smat)^(-1/q)*Smat
    if (p>0)
      beta <- beta[2:(p+1),]
    else 
      beta <- matrix(0,0,q)
      
    if (k>0) {
        # do the refining
        tmp <- resmultiGS(X,Y,beta,Cmat,k,b,cc,0,convTol)
        gammarw <- tmp$gammarw
        scalerw <- tmp$scalerw
        betarw <- tmp$betarw
        rdisrw <- tmp$rdis
    }
    
    else  # k = 0 means "no refining" 
       {  gammarw <- Cmat
        resrw <- res
        betarw <- beta
        places <- t(combn(1:n,2))
        term1vector <- as.matrix(res[places[,1],])
        term2vector <- as.matrix(res[places[,2],])
        diffrw <- term1vector-term2vector
       
        rdisrw <- sqrt(mahalanobis(diffrw, rep(0,q), gammarw))
        
        scalerw <- median(abs(rdisrw))/0.6745
    }
    
    if (i > 1) {   
        # if this isn't the first iteration....
        scaletest <- mean(rhobiweight(rdisrw/sworst,cc))
        if (scaletest < b) { 

            ss <- sort(bestscales, index.return=TRUE)
            ind <- ss$ix[bestr]
 
            bestscales[ind] <- scalemultiGS(rdisrw,b,cc,scalerw)
            bestbetas[,ind] <- vecop(betarw)
            bestgammas[,ind] <- vecop(gammarw)
            sworst <- max(bestscales)
        }
    }    
    else # if this is the first iteration, then this is the best beta...
       {
        bestscales[bestr] <- scalemultiGS(rdisrw,b,cc,scalerw)
        bestbetas[,bestr] <- vecop(betarw)
        bestgammas[,bestr] <- vecop(gammarw)

      }       
}
# do the complete refining step until convergence starting
# from the best subsampling candidate (possibly refined)
ibest <- which.min(bestscales)
superbestscale <- bestscales[ibest]
superbestbeta <- reconvec(bestbetas[,ibest],q)
superbestgamma <- reconvec(bestgammas[,ibest],q)


# magic number alert
for (i in bestr:1) { 
    tmp <- resmultiGS(X,Y,reconvec(bestbetas[,i],q), reconvec(bestgammas[,i],q),maxIt,b,cc,bestscales[i],convTol)
    if (tmp$scalerw < superbestscale)  {
        superbestscale <- tmp$scalerw
        superbestgamma <- tmp$gammarw
        superbestbeta <- tmp$betarw
    }
}
GSbeta <- superbestbeta
GScovariance <- superbestgamma*superbestscale^2
intercept <- IRLSlocation(Y-X%*%GSbeta,GScovariance,bdp=bdp,cc=cc)
GSbetatot <- rbind(intercept,GSbeta)

return(list( Beta = GSbetatot, Gamma = superbestgamma, scale = superbestscale, Sigma = GScovariance,b=b,c=cc))
}


#--------------------------------------------------------------------------------------

