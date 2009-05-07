plot.FRBpca <- function(x, which=1:3, pcs.loadings=1:min(5, length(x$eigval)), ...) {

    FRBres <- x
    currentAsk <- devAskNewPage(ask = NULL)

    if (is.element(1,which)) {
        plotFRBvars(FRBres, cumul=2)
        devAskNewPage(ask=TRUE)
    }
    if (is.element(2,which)) {
        plotFRBangles(FRBres)
        devAskNewPage(ask=TRUE)
    }
    if (is.element(3,which)) {
        plotFRBloadings(FRBres, pcs=pcs.loadings)
    }
    devAskNewPage(ask=currentAsk)

}
