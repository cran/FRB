summary.FRBmultireg <- function(object, confmethod = c("BCA","basic","both"), ...) {

confmethod <- match.arg(confmethod)

estBeta <- object$Beta
q <- ncol(estBeta)
p <- nrow(estBeta)
responses <- colnames(estBeta)
covariates <- rownames(estBeta)

if (!is.null(object$bootest))  {
  leftbr <- rep("(", p)
  rightbr <- rep(")", p)
  Betawstd <- as.data.frame(list(estBeta[,1], leftbr, object$SE[,1], rightbr))
  if (q>1) {
    for (i in 2:q) {
        Betawstd <- as.data.frame(list(Betawstd, estBeta[,i], leftbr, object$SE[,i], rightbr))
    }
  }
  colnames(Betawstd)[1+(0:(q-1))*4] <- responses
  colnames(Betawstd)[2+(0:(q-1))*4] <- " "
  colnames(Betawstd)[3+(0:(q-1))*4] <- " "
  colnames(Betawstd)[4+(0:(q-1))*4] <- " "
  
  if ((confmethod == "BCA") | (confmethod == "both")) {
    llist.bca <- list()
    for (i in 1:q) {
  
      limits <- cbind(object$CI.bca.lower[,i], object$CI.bca.upper[,i])
      rownames(limits) <- rownames(estBeta)
      colnames(limits) <- c("lower", "upper")
  
      significant <- sign(apply(limits,1,prod))
      starred <- rep("",p)
      starred[significant==1] <- "*"
      limits <- as.data.frame(list(limits, starred))
      colnames(limits)[3] <- " "
      
      llist.bca[[responses[i]]] <- limits
    }
  }
  if ((confmethod == "basic") | (confmethod == "both")) {
    llist.basic <- list()
    for (i in 1:q) {
  
      limits <- cbind(object$CI.basic.lower[,i], object$CI.basic.upper[,i])
      rownames(limits) <- rownames(estBeta)
      colnames(limits) <- c("lower", "upper")
  
      significant <- sign(apply(limits,1,prod))
      starred <- rep("",p)
      starred[significant==1] <- "*"
      limits <- as.data.frame(list(limits, starred))
      colnames(limits)[3] <- " "
  
      llist.basic[[responses[i]]] <- limits
    }
  }
  
  if (confmethod == "BCA")
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, limits.bca=llist.bca, method=object$method, conf=object$conf)
  else if (confmethod == "basic")
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, limits.basic=llist.basic, method=object$method, conf=object$conf)
  else
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, limits.bca=llist.bca, limits.basic=llist.basic, method=object$method, conf=object$conf)
}
else
    res <- list(responses=responses, covariates=covariates, Beta=estBeta, Sigma=object$Sigma, method=object$method)

class(res) <- "summary.FRBmultireg"

res

}

###############################################################################

print.summary.FRBmultireg <- function(x, digits=3, ...) {

if (x$method$est=="MM")
cat(paste("Multivariate regression based on multivariate MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n\n")
else
cat(paste("Multivariate regression based on multivariate ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n\n")

cat("Response variables: ", x$responses, "\n")
cat("Covariates: ", x$covariates, "\n\n")

if (!is.null(x$Betawstd)) {
  cat("Coefficient estimates (with bootstrap standard errors):\n")
  print(x$Betawstd, digits=digits)
} 
else {
  cat("Coefficient estimates):\n")
  print(x$Beta, digits=digits)
} 
if (!is.null(x$limits.bca)) {
  for (i in 1:length(x$limits.bca)) {
    cat("\n", x$conf*100, "% BCa confidence limits for response \"", x$responses[i], "\":\n", sep="")

    print(x$limits.bca[[x$responses[i]]], digits = digits)
  }
}

if (!is.null(x$limits.basic)) {
  for (i in 1:length(x$limits.basic)) {
    cat("\n", x$conf*100, "% \"basic bootstrap\" confidence limits for response \"", x$responses[i], "\":\n", sep="")

    print(x$limits.basic[[x$responses[i]]], digits = digits)
  }
}

cat("\n Error covariance matrix estimate:\n")
print(x$Sigma, digits=digits)

}

###############################################################################

print.FRBmultireg <- function(x, digits=3, ...) {


estBeta <- x$Beta
q <- ncol(estBeta)
p <- nrow(estBeta)

if (!is.null(x$bootest))  {
  leftbr <- rep("(", p)
  rightbr <- rep(")", p)
  Betawstd <- as.data.frame(list(estBeta[,1], leftbr, x$SE[,1], rightbr))
  if (q>1) {
    for (i in 2:q) {
        Betawstd <- as.data.frame(list(Betawstd, estBeta[,i], leftbr, x$SE[,i], rightbr))
    }
  }
  colnames(Betawstd)[1+(0:(q-1))*4] <- colnames(estBeta)
  colnames(Betawstd)[2+(0:(q-1))*4] <- " "
  colnames(Betawstd)[3+(0:(q-1))*4] <- " "
  colnames(Betawstd)[4+(0:(q-1))*4] <- " "
  
  if (x$method$est=="MM")
  cat(paste("Multivariate regression based on multivariate MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n\n")
  else
  cat(paste("Multivariate regression based on multivariate ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n\n")

  cat("Coefficient estimates (with bootstrap standard errors):\n")
  print(Betawstd, digits=digits)
}
else {
  if (x$method$est=="MM")
  cat(paste("Multivariate regression based on multivariate MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n\n")
  else
  cat(paste("Multivariate regression based on multivariate ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n\n")
  
  cat("Coefficient estimates:\n")
  print(estBeta, digits=digits)
}

}