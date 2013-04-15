\name{plot.FRBpca}
\alias{plot.FRBpca}
\alias{plotFRBvars}
\alias{plotFRBangles}
\alias{plotFRBloadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Method for Objects of class 'FRBpca' }
\description{
Plot functions for \code{FRBpca} objects: plots PC variances, PC angles and PC loadings, with bootstrap inference 
}
\usage{
\method{plot}{FRBpca}(x, which = 1:3, pcs.loadings = 1:min(5, length(x$eigval)),
confmethod = c("BCA","basic"), ...)


plotFRBvars(x, cumul = 2, confmethod = c("BCA","basic"), 
            npcs = min(10, length(x$eigval)))
plotFRBangles(x, pcs = 1:min(12,length(x$eigval)))
plotFRBloadings(x, confmethod = c("BCA","basic"), 
            pcs = 1:min(5, length(x$eigval)), nvars=min(10, length(x$eigval)))
}                
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an \R object of class \code{FRBpca}, typically created by \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}} }
  \item{which}{ integer number(s) between 1 and 3 to specify which plot is desired (1 = variances; 2 = angles; 3 = loadings)}
  \item{pcs.loadings}{ integer number(s) indicating for which of the PCs the loadings should be shown (in case the 
    \code{which} argument contains 2)}
  \item{cumul}{ integer between 0 and 2: 0 = screeplot, i.e. the variances of the PCs are shown; 1 = the cumulative variances (percentage) of the
      PCs are shown; 2 = (default) both plots are shown on the same page}
  \item{confmethod}{ which kind of bootstrap confidence intervals to be displayed: 
'BCA'= bias corrected and accelerated method, 'basic'= basic bootstrap method}
  \item{npcs}{ number of PCs to be included in screeplot/cumulative variances plot }
  \item{pcs}{ PCs to consider in plot; defaults to first 12 (maximally) for \code{plotFRBangles}; defaults to first 5 for
    \code{plotFRBloading} (each PC is on a separate page here) }
  \item{nvars}{ number of variables for which loadings should be shown in each PC; the loadings are shown in decreasing order in each PC}      
  \item{...}{ potentially more arguments } 
}
\details{
The generic \code{plot} function calls \code{plotFRBvars}, \code{plotFRBangles} and 
 \code{plotFRBloadings}, according to which
of these are respectively specified in argument \code{which}, and 
displays the plots on separate pages (the user is prompted for each new page).  The PCs for which the loadings should be plotted 
can be specified through the \code{pcs.loadings} argument. The other arguments are set to their default values by \code{plot}.

%Mainly for the purpose of outlier detection, \code{plotDiagPCA} shows the robust distances (i.e. mahalanobis distances based on the robust estimates of location and covariance).
%If \code{EIF=FALSE}, the distances are simply plotted against the index of the observations. 
%If \code{EIF=TRUE}, the diagnostic plot from Pison and Van Aelst (2002) is displayed.
%It plots the robust distance versus a measure of the overall empirical influence of the observation on the
%principal components. The empirical influences are obtained by using the influence function of
%the eigenvectors of the classical shape estimator at the normal model, and by
%substituting therein the robust estimates for the population parameters.
%The overall influence value is then defined by averaging the squared influence
%over all coefficients in the eigenvectors. 
%The vertical line on the plot is an indicative cutoff value, obtained through simulation. This last part takes 
%a few moments of computation time, which is the reason why the default in \code{plot} is \code{EIF=FALSE}.
%The default is reversed, however, when explicitly using the \code{plotDiagPCA} function. 

The solid curves displayed by \code{plotFRBvars} indicate the actual estimates of the variances (or percentages), while the dashed curves represent the confidence
limits as computed by \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}}.

\code{plotFRBangles} plots, for each PC,  histograms of the angles between the bootstrapped PC and the original PC estimate.
The angles are in radians, between 0 and pi/2. These limits are indicated by the red vertical lines. Angles close to zero correspond
to bootstrapped PCs closely aligned with the original PC, while an angle close to pi/2 means the bootstrapped PC is roughly
perpendicular to the original estimate (hence a large number of angles close to pi/2 implies high variability). 
If the number of PCs specified in \code{pcs} is very large (usually larger than the default settings), the histograms may not fit
on one page and a selection will be made (the user will be given a warning in that case).
 
In \code{plotFRBloadings}, the red dots represent the loadings, which are between -1 and 1. The square brackets indicate the confidence
limits as computed by \code{\link{FRBpcaS}} or \code{\link{FRBpcaMM}}. Only the loadings of the first \code{nvars} variables
are shown, where the variables were ordered according to the absolute value of the loading (i.e. only the \code{nvars} most important
variables for that particular PC are shown).
 
}
\references{ 
\itemize{
\item S. Van Aelst and G. Willems (2013). Fast and robust bootstrap for multivariate inference: The R package FRB. \emph{Journal of Statistical Software}, \bold{53}(3), 1--32.
        URL: \url{http://www.jstatsoft.org/v53/i03/}.
        }
}

\author{ Gert Willems and Ella Roelant }
\seealso{ \code{\link{FRBpcaS}}, \code{\link{FRBpcaMM}},  \code{\link{summary.FRBpca}} }
\examples{

data(ForgedBankNotes)

MMpcares <- FRBpcaMM(ForgedBankNotes, R=999, conf=0.95)
plot(MMpcares) 

# a closer look at the screeplot, specifying basic bootstrap intervals
plotFRBvars(MMpcares, cumul=0, confmethod="basic")

# plots the bootstrap angles for the first PC only
plotFRBangles(MMpcares, pcs=1)

# plots the loadings, with basic bootstrap intervals, for *all* the PCs 
plotFRBloadings(MMpcares, confmethod="basic", pcs=1:ncol(ForgedBankNotes))

}
