\name{MMest_loccov}
\alias{MMest_loccov}
\alias{Sest_loccov}
\alias{MMest_twosample}
\alias{Sest_twosample}
\title{S- and MM-Estimates of multivariate location and covariance matrix}
\description{
  Compute S- and MM-Estimates of multivariate location and covariance matrix 
}
\usage{
    MMest_loccov(Y, control=MMcontrol(...), ...)
    Sest_loccov(Y, bdp=.5, control=Scontrol(...), ...)
    MMest_twosample(X, groups, control=MMcontrol(...), ...)
    Sest_twosample(X, groups, bdp=0.5, control=Scontrol(...), ...)
}

\arguments{
  \item{Y}{input matrix or data frame}
  \item{X}{input matrix or data frame}
  \item{bdp}{breakdown point, defaults to 0.5}
  \item{groups}{grouping variable}
  \item{control}{a list with control parameters for tuning the S- or MM-estimate 
    and its computing algorithm, see\code{\link{Scontrol}} and  
    \code{\link{MMcontrol}}.}
  \item{...}{further arguments to be passed to \code{CovMMest()}}
}
\details{
    This functions are internal, wrappers around the functions \code{Sest()} \code{CovMMest()}.
}

\value{
Return lists with the following components:
  
  \item{Mu}{location}
  \item{Gamma}{shape}
  \item{scale}{scale=det^(1/(2*m))}
  \item{Sigma}{covariance matrix}
  \item{c1}{tuning parameter of the loss function for MM-estimation }
  \item{SMu}{location of the initial S-estimate}
  \item{SGamma}{shape of the initial S-estimate}
  \item{SSigma}{covariance matrix of the initial S-estimate}
  \item{b}{tuning parameters used in Tukey biweight loss function for S-estimation, as determined by bdp}
  \item{w}{scaled weights}
  \item{outflag}{outlier flags}
}
\examples{
\donttest{
    Y <- matrix(rnorm(50*5), ncol=5)
    (MMests <- MMest_loccov(Y)) 

    (Sests <- Sest_loccov(Y, bdp = 0.25)) 

    Y1 <- matrix(rnorm(50*5), ncol=5)
    Y2 <- matrix(rnorm(50*5), ncol=5)
    Ybig <- rbind(Y1,Y2)
    grp <- c(rep(1,50),rep(2,50))
    (MMests <- MMest_twosample(Ybig, grp))
}    
}
