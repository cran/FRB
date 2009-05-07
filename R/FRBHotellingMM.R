FRBhotellingMM <-function(Xdata,Ydata=NULL,mu0=0,R=999,conf=0.95,method=c("pool","HeFung"),control=MMcontrol(...), ...)
{
# performs robust Hotelling test based on multivariate MM estimates 
# with fast and robust bootstrap
#
# calls: MMest_loccov(), MMboot_loccov(), twosampleMM(), MMboottwosample()
#
# Input
# Xdata: (n x p) data set
# mu0: mean under null hypothesis in case of one sample
# R: number of bootstrap samples
# method: "pool" uses the pooled covariance matrix  
#          "HeFung" uses the covariance matrix of He and Fung
#          to estimate the common covariance matrix 
# Output
# p-value from the robust Hotelling test
#-------------------------------------------------------------------------
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

#--------------------------------------------------------------------------
hottest1MM<-function(Xdata,mu0,R,conf,control=MMcontrol())
{
#Robust one-sample Hotelling test based on MM-estimator

    mu0<- as.vector(mu0)
    n <- nrow(Xdata)
    p <- ncol(Xdata)
    dimens <- p+p*p
    

    MMests <- MMest_loccov(Xdata,control=control) 
    MMmu <- MMests$Mu
    MMSigma <- MMests$Sigma
    MMGamma <- MMests$Gamma
    Smu <- MMests$SMu
    SSigma <- MMests$SSigma
    teststatMM<-n*(MMmu -mu0) %*%solve(MMSigma)%*%t(MMmu-mu0)
    testvectMM<-c()
    bootresMM<-MMboot_loccov(Xdata,R=R,ests=MMests)
    bmatrixSenmm <- bootresMM$centered

for (r in 1:R) {
        estimates<-bmatrixSenmm[,r]
        estmubminmuhat<-estimates[1:p]
	  estmubminmuhat<-as.matrix(estmubminmuhat)	
        estgammabmingammahat<-estimates[(p+1):dimens]
        estgammab<-estgammabmingammahat+vecop(MMGamma)
        estgammab<-reconvec(estgammab,p)
        estsigmabminsigmahat<-estimates[(dimens+1):(dimens+p*p)]
        estsigmab<-estsigmabminsigmahat+vecop(SSigma)
        estsigmab<-reconvec(estsigmab,p)
        if ((det(estsigmab) >0) && (det(estgammab) >0)) {
            estsigmammb<-det(estsigmab)^(1/p)*estgammab
            testvectMMwaarde<-n*t(estmubminmuhat)%*%solve(estsigmammb)%*%estmubminmuhat
            testvectMM<-c(testvectMM,testvectMMwaarde)
        }    
}    
pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok <- length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
	   conf.int[1,i] <- MMmu[i]-sqrt(1/n*quantval*MMSigma[i,i])
     conf.int[2,i] <- MMmu[i]+sqrt(1/n*quantval*MMSigma[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(MMmu) <- ("   Estimate")


return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc=MMmu,cov=MMSigma,confint=conf.int))

}

#-----------------------------------------------------------------------------

hottest2MMpool <-function(Xdata1,Xdata2,R,conf,control=MMcontrol())
{
#Robust two-sample Hotelling test based on MM-estimator and pooled covariance matrix


n1 <- nrow(Xdata1)
p <- ncol(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
dimens<-p+p*p
MMests1 <- MMest_loccov(Xdata1,control=control)
MMests2 <- MMest_loccov(Xdata2,control=control)
MMmu1 <- MMests1$Mu
MMSigma1 <- MMests1$Sigma
MMGamma1 <- MMests1$Gamma
Smu1<- MMests1$SMu
SSigma1 <- MMests1$SSigma
MMmu2 <- MMests2$Mu
MMSigma2 <- MMests2$Sigma
Smu2<-MMests2$SMu
SSigma2<-MMests2$SSigma
MMGamma2 <- MMests2$Gamma

MMSigmap<-((n1-1)*MMSigma1+(n2-1)*MMSigma2)/(n1+n2-2)
SSigmap<-((n1-1)*SSigma1+(n2-1)*SSigma2)/(n1+n2-2)
teststatMM<-((n1*n2)/(n1+n2))*(MMmu1-MMmu2)%*%solve(MMSigmap)%*%t(MMmu1-MMmu2)


bootresMM1<-MMboot_loccov(Xdata1,R=R,ests=MMests1)
bootresMM2<-MMboot_loccov(Xdata2,R=R,ests=MMests2)
bmatrixSenmm1 <- bootresMM1$centered
bmatrixSenmm2 <- bootresMM2$centered

testvectMM<-c()

    
for (r in 1:R) {
    estimates1<-bmatrixSenmm1[,r]
    estmubminmuhat1<-estimates1[1:p]
    estmubminmuhat1 <- as.matrix(estmubminmuhat1)
    estgammabmingammahat1<-estimates1[(p+1):dimens]
    estgammab1<-estgammabmingammahat1+vecop(MMGamma1)
    estgammab1<-reconvec(estgammab1,p)
    estsigmabminsigmahat1<-estimates1[(dimens+1):(dimens+p*p)]
    estsigmab1<-estsigmabminsigmahat1+vecop(SSigma1)
    estsigmab1<-reconvec(estsigmab1,p)

    estimates2<-bmatrixSenmm2[,r]
    estmubminmuhat2<-estimates2[1:p]
    estmubminmuhat2 <- as.matrix(estmubminmuhat2)
    estgammabmingammahat2<-estimates2[(p+1):dimens]
    estgammab2<-estgammabmingammahat2+vecop(MMGamma2)
    estgammab2<-reconvec(estgammab2,p)
    estsigmabminsigmahat2<-estimates2[(dimens+1):(dimens+p*p)]
    estsigmab2<-estsigmabminsigmahat2+vecop(SSigma2)
    estsigmab2<-reconvec(estsigmab2,p)

   

    if ((det(estsigmab1) >0) && (det(estsigmab2)>0) && (det(estgammab1) >0) && (det(estgammab2)>0)) {
        estsigmammb1<-det(estsigmab1)^(1/p)*estgammab1
        estsigmammb2<-det(estsigmab2)^(1/p)*estgammab2
        estsigmaMMP<-((n1-1)*estsigmammb1+(n2-1)*estsigmammb2)/(n1+n2-2)
        testvectMMwaarde<-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaMMP)%*%(estmubminmuhat1-estmubminmuhat2)

        testvectMM=c(testvectMM,testvectMMwaarde)
    }
}

pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok<-length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
	   conf.int[1,i] <- MMmu1[i]-MMmu2[i]-sqrt(1/n*quantval*MMSigmap[i,i])
	   conf.int[2,i] <- MMmu1[i]-MMmu2[i]+sqrt(1/n*quantval*MMSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(MMmu1) <- rownames(MMmu2) <- ("   Estimate")



return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc1=MMmu1,loc2=MMmu2,cov=MMSigmap,confint=conf.int))


}
#----------------------------------------------------------------------------

hottest2MMHe<-function(Xdata1,Xdata2,R,conf,control=MMcontrol()){
#Robust two-sample Hotelling test based on MM-estimator of He and Fung 

n1 <- nrow(Xdata1)
n2 <- nrow(Xdata2)
n <- n1*n2/(n1+n2)
p <- ncol(Xdata1)
dimens <- 2*p+p*p
Xdata <- rbind(Xdata1,Xdata2)
groups <- c(rep(1,n1),rep(2,n2))
estsMM <- twosampleMM(Xdata,groups,control=control)

MMmu1 <- estsMM$Mu1
Smu1 <- estsMM$SMu1
MMmu2 <- estsMM$Mu2
Smu2 <- estsMM$SMu2
SSigmap <- estsMM$SSigma
MMGammap <- estsMM$Gamma
MMSigmap <- estsMM$Sigma

teststatMM <-((n1*n2)/(n1+n2))*(MMmu1-MMmu2)%*%solve(MMSigmap)%*%t(MMmu1-MMmu2)

bootresMM <- MMboottwosample(Xdata, groups=groups, R, ests=estsMM)
bmatrixSenmm <- bootresMM$centered
testvectMM <- c()

for (r in 1:R) {
    estimates <- bmatrixSenmm[,r]
    estmubminmuhat1 <- estimates[1:p]
    estmubminmuhat2 <- estimates[(p+1):(2*p)]
    estmubminmuhat1 <- as.matrix(estmubminmuhat1)
    estmubminmuhat2 <- as.matrix(estmubminmuhat2)
    estgammabmingammahat <- estimates[(2*p+1):(2*p+p*p)]
    estgammaMMP <- estgammabmingammahat+vecop(MMGammap)
    estgammaMMP <- reconvec(estgammaMMP,p)
    estsigmabminsigmahatSP <- estimates[(2*p+p*p+1):(2*p+2*p*p)]
    estsigmabSP <- estsigmabminsigmahatSP+vecop(SSigmap)
    estsigmabSP <- reconvec(estsigmabSP,p)
    if ((det(estsigmabSP) >0) & (det(estgammaMMP) >0)) {
        estsigmaMMP <- det(estsigmabSP)^(1/p)*estgammaMMP
        testvectMMwaarde<-((n1*n2)/(n1+n2))*t(estmubminmuhat1-estmubminmuhat2)%*%solve(estsigmaMMP)%*%(estmubminmuhat1-estmubminmuhat2)
        testvectMM=c(testvectMM,testvectMMwaarde)
    }
}
    

pvalue<-mean(testvectMM>=as.numeric(teststatMM))
Rok<-length(testvectMM)
testquantile <- floor((1 - (1 - conf)/2) * Rok)
quantval <- testvectMM[testquantile]
conf.int <- matrix(1:(2*p), ncol=p)
for (i in 1:p) {
    conf.int[1,i] <- MMmu1[i]-MMmu2[i]-sqrt(1/n*quantval*MMSigmap[i,i])
		conf.int[2,i] <- MMmu1[i]-MMmu2[i]+sqrt(1/n*quantval*MMSigmap[i,i])
}
dimnames(conf.int) <- list(c("Lower bound","Upper bound"),dimnames(Xdata)[[2]]) 
rownames(MMmu1) <- rownames(MMmu2) <- ("   Estimate")


return(list(teststat=teststatMM,testvect=testvectMM,pvalue=pvalue,loc1=MMmu1,loc2=MMmu2,cov=MMSigmap,confint=conf.int))


}

#----------------------------------------------------------------------------
#-                         main function                                     -
#----------------------------------------------------------------------------

method <- match.arg(method)

p <- ncol(Xdata)
n <- nrow(Xdata)

Xdatam <- as.matrix(Xdata)
if (is.null(colnames(Xdatam)))
    colnames(Xdatam) <- paste("V",1:p,sep="")

if (is.null(Ydata)==TRUE)
 {nrsamples <- 1}
else
 {Ydatam <- as.matrix(Ydata)
  if (is.null(colnames(Ydatam)))
    colnames(Ydatam) <- paste("V",1:p,sep="")
  nrsamples <- 2}

if (length(mu0)==1)
      {mu0=rep(mu0,p)}

if (nrsamples==1)
   {
   meth=paste("One sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")", sep="")
   res=hottest1MM(Xdatam,mu0=mu0,R=R,conf=conf,control=control)
   z <- list(pvalue=res$pvalue,teststat=res$teststat,teststat.boot=res$testvect,Mu=res$loc,Sigma=res$cov,
            CI=res$confint,Mu0=mu0,conf=conf,data=substitute(Xdata),meth=meth)
}
else 
    {
    if (method == "pool") {
      meth=paste("Two sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")\n",
          "(common covariance estimated by pooled covariance matrix)", sep="")
      res=hottest2MMpool(Xdatam,Ydatam,R=R,conf=conf,control=control)
    }
    else { 
      meth=paste("Two sample Hotelling test based on multivariate MM-estimates (bdp = ", control$bdp,", eff = ", control$eff, ")\n",
          "(common covariance estimated by He and Fung method)", sep="")
      res=hottest2MMHe(Xdatam,Ydatam,R=R,conf=conf,control=control)  
    }
    z <- list(pvalue=res$pvalue,teststat=res$teststat,teststat.boot=res$testvect,Mu1=res$loc1,Mu2=res$loc2,Sigma=res$cov,
            CI=res$confint,conf=conf,data=c(substitute(Xdata),substitute(Ydata)),meth=meth)
}

class(z) <- "FRBhot"

return(z)

}  
#-----------------------------------------------------------------------------------------









