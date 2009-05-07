\name{ForgedBankNotes}
\alias{ForgedBankNotes}
\docType{data}
\title{ Swiss (forged) bank notes data}
\description{
  Six measurements made on 100 forged Swiss bank notes. 
}
\usage{data(ForgedBankNotes)}
\format{
  The data frame contains the following columns:
  
% chr "ForgedBankNotes"
\describe{
\item{\code{Length}}{metric length of the bill}
\item{\code{Left}}{height of the bill, measured on the left}
\item{\code{Right}}{height of the bill, measured on the right}
\item{\code{Bottom}}{distance of inner frame to the lower border}
\item{\code{Top}}{distance of inner frame to the upper border}
\item{\code{Diagonal}}{length of the diagonal}
}
}

\details{
  The original data set in Flury and Riedwyl (1988) additionally contained 100 genuine bank notes, but these are not included here. 
}
\source{
 B. Flury and H. Riedwyl (1988). \emph{Multivariate Statistics: A practical approach}. London: Chapman & Hall. 
 }
\references{
M. Salibian-Barrera, S. Van Aelst and G. Willems (2006). PCA based on multivariate MM-estimators
with fast and robust bootstrap. \emph{Journal of the American Statistical Association}, \bold{101}, 1198-1211.
}
\examples{
data(ForgedBankNotes)
pairs(ForgedBankNotes)
%## maybe str(ForgedBankNotes) ; plot(ForgedBankNotes) ...
}
\keyword{datasets}
