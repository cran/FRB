\name{MMest_loccov}
\alias{MMest_loccov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MM-Estimates of Location and Covariance }
\description{
  Computes MM-estimates of multivariate location and covariance, using an initial S-estimate
}
\usage{
MMest_loccov(Y, control=MMcontrol(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ matrix or data frame }
  \item{control}{a list with control parameters for tuning the MM-estimate and its computing algorithm, 
  see \code{\link{MMcontrol}}().}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}}
}
\details{
  This function is called by \code{\link{FRBpcaMM}} and \code{\link{FRBhotellingMM}}. 
  
  The MM-estimates are defined by first computing S-estimates of location and covariance, then fixing the scale component and 
  re-estimating the location and shape by more efficient M-estimates (see Tatsuoka and Tyler (2000)). Tukey's biweight is used for 
  the loss functions. By default the first loss function (in the S-estimates) is tuned in order to obtain 50\% breakdown point. 
  The default tuning of the second loss function (M-estimates) ensures 95\% efficiency at the normal model. 
  This tuning can be changed via argument \code{control} if desired. 
  When interested in location estimates, control parameter \code{shapeEff} should be \code{FALSE} (the default), 
  in which case the particular efficiency is that of the location estimate. When interest lies in the covariance or shape part, e.g. in PCA analysis, 
  it makes sense to set \code{shapeEff=TRUE}, in which case the shape efficiency is considered instead.
      
  The computation of the S-estimates is performed by a call to \code{\link{Sest_loccov}}, which uses the fast-S algorithm. 
  See \code{\link{MMcontrol}}() to see or change the tuning parameters for this algorithm.  The M-estimate part is computed
  through iteratively reweighted least squares (RWLS).
  
  Apart from the MM-location estimate \code{Mu}, the function returns both the MM-covariance \code{Sigma} and MM-shape estimate \code{Gamma} (which has determinant equal to 1). 
  Additionally, the S-estimates are returned as well (their Gaussian efficiency is usually lower than the MM-estimates but they may 
  have a lower bias). 
}
\value{
A list containing:
  \item{Mu }{MM-estimate of location}
  \item{Sigma }{MM-estimate of covariance}
  \item{Gamma }{MM-estimate of shape}
  \item{SMu }{S-estimate of location}
  \item{SSigma }{S-estimate of covariance}
  \item{SGamma }{S-estimate of shape}
  \item{scale }{S-estimate of scale (univariate)}
  \item{c0,b,c1}{tuning parameters of the loss functions (depend on control parameters \code{bdp} and \code{eff})}
  \item{w}{implicit weights corresponding to the MM-estimates (i.e. final weights in the RWLS procedure)}
  \item{outFlag}{outlier flags: 1 if the robust distance of the observation exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of \code{Y}; 0 otherwise}
}
\references{
\itemize{ 
\item K.S. Tatsuoka and D.E. Tyler (2000). The uniqueness of S and M-functionals under non-elliptical distributions.
\emph{The Annals of Statistics}, \bold{28}, 1219-1243. }
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{Sest_loccov}}, \code{\link{FRBpcaMM}}, \code{\link{FRBhotellingMM}},  
\code{\link{MMboot_loccov}}, \code{\link{MMcontrol}} }
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
Y <- matrix(rnorm(50*5), ncol=5)
MMests <- MMest_loccov(Y) 
# MM-estimate of location:
MMests$Mu
# MM-estimate of covariance:
MMests$Sigma
# initial S-estimate of location:
MMests$SMu

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
