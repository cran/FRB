\name{plot.FRBmultireg}
\alias{plot.FRBmultireg}
\alias{plotDiag}
\alias{plotFRBconf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Method for Objects of class 'FRBmultireg' }
\description{
  Plot functions for objects of class \code{FRBmultireg}. Two types of plots can be produced: (1) a diagnostic plot, presenting
  (multivariate) residuals versus robust distances in the covariates space; and (2) histograms for the bootstrap estimates for all 
  (or a selection) of the regression coefficients, based on Fast and Robust Bootstrap and with visualization of bootstrap confidence limits.
}
\usage{
\method{plot}{FRBmultireg}(x, which=1:2, ...)

plotDiag(x, onepage = TRUE, ...)
plotFRBconf(x, expl, resp, confmethod = c("BCA","basic"), onepage = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBmultireg}, typically created by  \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}} or \code{\link{FRBmultiregGS}} }
  \item{which}{ integer number(s) between 1 and 2 to specify which plot is desired  }
  \item{expl}{ optional; vector specifying the explanatory variables to be shown (either by index or by variable name) }
  \item{resp}{ optional; vector specifying the response variables to be shown (either by index or by variable name) }
  \item{confmethod}{ which kind of bootstrap confidence intervals to be displayed: 'BCA'= bias corrected and accelerated method, 
                  'basic'= basic bootstrap method}
  \item{onepage}{ logical: if TRUE, all requested histograms are plotted on one page; if FALSE, separate pages are used for
                  each response variable } 
  \item{\dots}{ potentially more arguments to be passed }
}
\details{
  The generic \code{plot} function calls \code{plotDiag} and \code{plotFRBconf}, according to which
of these are respectively specified in argument \code{which}, and displays the plots on separate pages (the user is prompted for each new page).

The diagnostic plot, obtained through \code{plotDiag}, shows the residual distances based on the robust coefficient and covariance
estimates in \code{x}, versus the distances within the covariates space. The latter are based on robust estimates of location and scatter for the
data matrix \code{x$X} (without intercept). Computing these robust estimates may take an appreciable amount of time. The estimator used
corresponds to the one which was used in obtaining \code{x} (with the same breakdown point, for example, and the same control parameters).
On the vertical axis a cutoff line is drawn at the square root of the .975 quantile of the chi-squared distribution with degrees of 
freedom equal to the number of response variables. On the horizontal axis the same quantile is drawn but now with degrees of freedom
equal to the number of covariates (not including intercept).
Those points to the right of the cutoff can be viewed as high-leverage points. These can be classified into so-called 
'bad' or 'good' leverage points depending on whether they are above or below the cutoff. Points above the cutoff but to the
left of the vertical cutoff are sometimes called vertical outliers.
See also Van Aelst and Willems (2005) for example.
 
With \eqn{p} and \eqn{q} the number of explanatory resp. response variables specified, \code{plotFRBconf} by default 
  (i.e. if \code{onepage=TRUE})
  plots a \eqn{p} by \eqn{q} matrix of histograms, showing the bootstrap recalculations of the corresponding entry in the regression 
  coefficient matrix \code{Beta} as provided in \code{x}. The original estimates for the coefficients are indicated by dotted lines,
  while the solid lines are the bootstrap confidence limits. In case the interval does not contain zero, the plot title is
  printed in red, indicating significance.
  
%  However, if \eqn{q>5} or \eqn{p>4} (and \code{expl} and \code{resp} are not provided), then a subset of coefficients to plot is selected
%  since otherwise the number of plots may be too large to fit on one page, which may result in an error. 
%  If \code{expl} or \code{resp} are provided, no subset is selected. In that case, and in (default) mode 
%  \code{onepage=TRUE}, select no more than 5 variables each, to avoid errors. Or use \code{onepage=FALSE}.  
  However, if \eqn{p} and/or \eqn{q} are large, the histograms may not fit on the page and an attempt to do it may
  result in an error. Therefore, the function first tries whether it fits (the outcome is platform-dependent), 
  and if not it reduces \eqn{p} and/or \eqn{q} until all plots do fit on the page. Hence, only a selection may be shown and the
  user is given a warning in that case.   

  If \code{onepage=FALSE}, separate pages are used for each response variable and the user is prompted for page change.
  In case the number (\eqn{p}) of explanatory variables is very large, the function again may show only a selection. 
  
  The arguments for \code{plotFRBconf} may also be passed throug the generic \code{plot} function.
}
%\value{
%  Produces one or more arrays of histograms.
%}
\references{ 
\itemize{
\item  S. Van Aelst and G. Willems (2005). Multivariate regression S-estimators for robust estimation and     libra
inference. \emph{Statistica Sinica}, \bold{15}, 981-1001.
}
}
\author{ Gert Willems and Ella Roelant }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{FRBmultiregS}}, \code{\link{FRBmultiregMM}}, \code{\link{FRBmultiregGS}}, \code{\link{summary.FRBmultireg}}  }
\examples{
data(schooldata)
school.x <- data.matrix(schooldata[,1:5])
school.y <- data.matrix(schooldata[,6:8])

Sres <- FRBmultiregS(school.x, school.y, R=999, bdp = 0.25, conf = 0.99)

plot(Sres)
# both types of plots are produced
# computing the robust distances takes some time

#if only interested in the confidence intervals, this is faster:
plot(Sres, which=2)
# or, equivalently,
plotFRBconf(Sres)
#the significant coefficients can be spotted as those with red plot title

#the plotFRBconf command above selected a subset, since otherwise an error may occur; 
#as may happen when you explicitely ask for all coefficients to be plotted on one page:
\dontrun{
plot(Sres, which=2, expl=1:6, resp=1:3)
}

#use separate pages for each response in case of many covariates: 
plot(Sres, which=2, onepage=FALSE)

#perhaps specify some specific variables of interest:
plot(Sres, which=2, expl=c("education", "occupation"), resp=c("selfesteem","reading"))
#or (the same):
plot(Sres, which=2, expl=2:3, resp=c(3,1))


}
