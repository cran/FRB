plot.FRBmultireg <- function(x, which=1:2, ...) {

currentAsk <- devAskNewPage(ask = NULL)

if (is.element(1,which)) {
    plotDiag(x,...)
    devAskNewPage(ask=TRUE)
}
if (is.element(2,which)) {
    plotFRBconf(x,...)
    devAskNewPage(ask=TRUE)
}
devAskNewPage(ask = currentAsk) 

}

