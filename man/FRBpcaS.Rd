\name{FRBpcaS}
\alias{FRBpcaS}
\alias{FRBpcaS.default}
\alias{FRBpcaS.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ PCA based on Multivariate S-estimators with Fast and Robust Bootstrap }
\description{
  Performs principal components analysis based on the robust S-estimate of the shape matrix. Additionally uses the Fast and Robust Bootstrap
  method to compute inference measures such as standard errors and confidence intervals. 
}
\usage{
\method{FRBpcaS}{formula}(formula, data=NULL, ...)

\method{FRBpcaS}{default}(Y, R = 999, bdp = 0.5, conf = 0.95, control=Scontrol(...),
na.action=na.omit, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}}; a symbolic description of the model to be fit.}
  \item{data}{ data frame from which variables specified in formula are to be taken.}
  \item{Y}{ matrix or data frame. }
  \item{R}{ number of bootstrap samples. Default is \code{R=999}. }
  \item{bdp}{ required breakdown point for the S-estimates. Should have \eqn{0 < } \code{bdp} \eqn{\le 0.5}, the default is 0.5. }
  \item{conf}{level of the bootstrap confidence intervals. Default is \code{conf=0.95}. }
  \item{control}{a list with control parameters for tuning the computing algorithm, see \code{\link{Scontrol}}().}
  \item{na.action}{ a function which indicates what should happen when the data contain NAs. Defaults to \code{\link{na.omit}}.}
  \item{...}{allows for specifying control parameters directly instead of via \code{control}.}
}
\details{                                     
  Multivariate S-estimates were introduced by Davies (1987) and can be highly robust while enjoying a reasonable Gaussian efficiency.
  The loss function used here is Tukey's biweight. It will be tuned in order to achieve the required breakdown point \code{bdp} 
(any value between 0 and 0.5). 
  The MM-estimates are computed by a call to the implementation of the fast-S algorithm (Salibian-Barrera and Yohai 2006) in the \pkg{rrcov} package of Todorov and Filzmoser (2009). \code{\link{Scontrol}} provides some adjustable tuning parameters regarding the algorithm. The result of this call is also returned as the value \code{est}.
  
  PCA is performed by computing the eigenvalues (\code{eigval}) and eigenvectors (\code{eigvec}) of the S-estimate of shape, which is a rescaled version of the
  S-estimate of covariance (rescaled to have determinant equal to 1). With \code{pvar} the function also provides the estimates for the percentage of 
  variance explained by the first \eqn{k} principal components, which are simply the cumulative proportions of the eigenvalues sum. 
  Here, \eqn{k} ranges from 1 to \eqn{p-1} (with \eqn{p} the number of variables in \code{Y}). 
  The eigenvectors are always given in the order of descending eigenvalues.
  
  The Fast and Robust Bootstrap (Salibian-Barrera and Zamar 2002) is used to calculate standard errors, and also so-called
  basic bootstrap confidence intervals and bias corrected and accelerated (BCa) confidence intervals (Davison and Hinkley 1997, p.194 and p.204 respectively) corresponding
  to the estimates \code{eigval}, \code{eigvec} and \code{pvar}. 
  The bootstrap is also used to estimate the average angles between true and estimated eigenvectors, returned as \code{avgangle}.
  See Salibian-Barrera, Van Aelst and Willems (2006).
  The fast and robust bootstrap computations for the S-estimates are performed by \code{\link{Sboot_loccov}}() and its raw result can be found in \code{bootest}.
  The actual bootstrap values of the PCA-related quantities can be found in \code{eigval.boot}, \code{eigvec.boot} and \code{pvar.boot}, where each column 
  represents a bootstrap sample. For \code{eigvec.boot}, the eigenvectors are stacked on top of each other and the same goes for 
  \code{eigvec.CI.bca} and \code{eigvec.CI.basic} which hold the confidence limits. 
  
  The two columns in the confidence limits always respectively represent the lower and upper limits.
  For the percentage of variance the function also provides one-sided confidence intervals ([-infty upper]), which can be
  used to test the hypothesis that the true percentage at least equals a certain value.  

  Bootstrap samples are discarded if the fast and robust covariance estimate is not positive definite, such that the actual number 
  of recalculations used can be lower than \code{R}. This actual number equals \code{R} - \code{failedsamples}. However, if 
  more than \code{0.75R} of the bootstrap shape estimates is non-positive definite, the failed bootstrap samples are recovered 
  by applying the \code{make.positive.definite} function (from package \code{corpcor}). If this also fails, the corresponding bootstrap sample
  is discarded after all, but such situation should be rare.
  This recovery may have an impact on the confidence limits and standard errors of especially the smallest eigenvalues in 
  \code{eigval} and \code{pvar}.
  }
\value{
  An object of class \code{FRBpca}, which contains the following components:
\item{shape}{(p x p) S-estimate of the shape matrix of \code{Y}}
\item{eigval}{(p x 1) eigenvalues of S shape }
\item{eigvec}{(p x p) eigenvectors of S-shape }
\item{pvar}{(p-1 x 1) percentages of variance for S eigenvalues }
\item{eigval.boot}{(p x R) eigenvalues of S shape }
\item{eigvec.boot}{(p*p x R) eigenvectors of S-shape (vectorized) }
\item{pvar.boot}{(p-1 x R) percentages of variance for S eigenvalues }
\item{eigval.SE}{(p x 1) bootstrap standard error for S eigenvalues }
\item{eigvec.SE}{(p x p) bootstrap standard error for S eigenvectors }
\item{pvar.SE}{(p-1 x 1) bootstrap standard error for percentage of variance for S eigenvalues }
\item{angles}{(p x R) angles between bootstrap eigenvectors and original S eigenvectors (in radians; in [0 pi/2])}
\item{avgangle}{(p x 1) average angles between bootstrap eigenvectors and original S eigenvectors (in radians; in [0 pi/2])}
\item{eigval.CI.bca}{(p x 2) BCa intervals for S eigenvalues }
\item{eigvec.CI.bca}{(p*p x 2) BCa intervals for S eigenvectors (vectorized)}
\item{pvar.CI.bca}{(p-1 x 2) BCa intervals for percentage of variance for S-eigenvalues}
\item{pvar.CIone.bca}{(p-1 x 1) one-sided BCa intervals for percentage of variance for S-eigenvalues ([-infty upper]) } 
\item{eigval.CI.basic}{(p x 2) basic bootstrap intervals for S eigenvalues }
\item{eigvec.CI.basic}{(p*p x 2) basic bootstrap intervals for S eigenvectors (vectorized)}
\item{pvar.CI.basic}{(p-1 x 2) basic bootstrap intervals for percentage of variance for S-eigenvalues }
\item{pvar.CIone.basic}{(p-1 x 1) one-sided basic bootstrap intervals for percentage of variance for S-eigenvalues ([-infty upper]) }
\item{est}{list containing the S-estimates of location and scatter }
\item{bootest}{(list) result of \code{\link{Sboot_loccov}}() }
\item{failedsamples}{number of bootstrap samples with non-positive definiteness of shape }
\item{conf}{ a copy of the \code{conf} argument}
\item{method}{ a character string giving the robust PCA method that was used}
\item{w}{implicit weights corresponding to the S-estimates (i.e. final weights in the RWLS procedure at the end of the fast-S algorithm)}
\item{outFlag}{outlier flags: 1 if the robust distance of the observation exceeds the .975 quantile of (the square root of)
  the chi-square distribution with degrees of freedom equal to the dimension of \code{Y}; 0 otherwise}
\item{Y}{copy of the data argument as a matrix}
}
\references{
\itemize{
\item P.L. Davies (1987) Asymptotic behavior of S-estimates of multivariate location parameters and dispersion
matrices. \emph{The Annals of Statistics}, \bold{15}, 1269-1292.
\item A.C. Davison and D.V. Hinkley (1997) \emph{Bootstrap Methods and their Application}. Cambridge Series in
Statistical and Probabilistic Mathematics. Cambridge: Cambridge University Press.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2006) PCA based on multivariate MM-estimators
with fast and robust bootstrap. \emph{Journal of the American Statistical Association}, \bold{101}, 1198-1211.
\item M. Salibian-Barrera, S. Van Aelst and G. Willems (2008) Fast and robust 
bootstrap. \emph{Statistical Methods and Applications}, \bold{17}, 41-71. 
\item M. Salibian-Barrera, R.H. Zamar (2002) Bootstrapping robust estimates of 
regression. \emph{The Annals of Statistics}, \bold{30}, 556-582. 
\item V. Todorov and P. Filzmoser (2009) An object-oriented framework for robust multivariate analysis. \emph{Journal of Statistical Software}, \bold{32}, 1--47. URL \url{http://www.jstatsoft.org/v32/i03/}.
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}
\author{ Gert Willems, Stefan Van Aelst and Ella Roelant }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{plot.FRBpca}}, \code{\link{summary.FRBpca}}, \code{\link{print.FRBpca}}, \code{\link{FRBpcaMM}}, 
\code{\link{Sboot_loccov}}, \code{\link{Scontrol}} }
\examples{
 
data(ForgedBankNotes)

Spcares <- FRBpcaS(ForgedBankNotes, R=999, bdp=0.25, conf=0.95)
# or using the formula interface
\dontrun{Spcares <- FRBpcaMM(~.,data=ForgedBankNotes, R=999, conf=0.95)}

# the simple print method shows the standard deviations with confidence limits:
Spcares

# the summary functions shows a lot more (see help(summary.FRBpca)):
summary(Spcares)

# ask for the eigenvalues:
Spcares$eigval

# or, in more pretty format, with confidence limits:
summary(Spcares)$eigvals

# note that the standard deviations of the print-output can also be asked for by:
sqrt( summary(Spcares)$eigvals )

# the eigenvectors and their standard errors:
Spcares$eigvec   # or prettier: summary(MMpcares)$eigvecs
Spcares$eigvec.SE
 
# take a look at the bootstrap distribution of the first eigenvalue
hist(Spcares$eigval.boot[1,])
# that bootstrap distribution is used to compute confidence limits as depicted 
# by the screeplot function:
plotFRBvars(Spcares, cumul=0)

# all plots for the FRB-PCA result:
plot(Spcares)
}

\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ robust }
