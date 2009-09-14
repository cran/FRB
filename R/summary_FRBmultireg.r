summary.FRBmultireg <- function(object, confmethod = c("BCA","basic","both"), digits=3, ...) {

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
  
  separator <- rep("|", p)
  sigcodes.bca <- matrix(rep("",p*q), ncol=q)
  sigcodes.bca[object$p.bca<0.1] <- "."  
  sigcodes.bca[object$p.bca<0.05] <- "*"  
  sigcodes.bca[object$p.bca<0.01] <- "**"  
  sigcodes.bca[object$p.bca<0.001] <- "***"  
  sigcodes.basic <- matrix(rep("",p*q), ncol=q)
  sigcodes.basic[object$p.basic<0.1] <- "."  
  sigcodes.basic[object$p.basic<0.05] <- "*"  
  sigcodes.basic[object$p.basic<0.01] <- "**"  
  sigcodes.basic[object$p.basic<0.001] <- "***"  
  if ((confmethod == "BCA") | (confmethod == "both")) {
    llist.bca <- list()
    for (i in 1:q) {
  
#      limits <- cbind(object$CI.bca.lower[,i], object$CI.bca.upper[,i])
#      rownames(limits) <- rownames(estBeta)
#      colnames(limits) <- c("lower", "upper")
      BetaTable <- as.data.frame(list(estBeta[,i], separator, object$SE[,i], separator, object$CI.bca.lower[,i], object$CI.bca.upper[,i], separator, object$p.bca[,i], sigcodes.bca[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "lower", "upper", " ", "p-value", " ")
  
#      significant <- sign(apply(limits,1,prod))
#      starred <- rep("",p)
#      starred[significant==1] <- "*"
#      limits <- as.data.frame(list(limits, starred))
#      colnames(limits)[3] <- " "
      
      llist.bca[[responses[i]]] <- BetaTable
    }
  }
  if ((confmethod == "basic") | (confmethod == "both")) {
    llist.basic <- list()
    for (i in 1:q) {
  
#      limits <- cbind(object$CI.basic.lower[,i], object$CI.basic.upper[,i])
#      rownames(limits) <- rownames(estBeta)
#      colnames(limits) <- c("lower", "upper")

      BetaTable <- as.data.frame(list(estBeta[,i], separator, object$SE[,i], separator, object$CI.basic.lower[,i], object$CI.basic.upper[,i], separator, object$p.basic[,i], sigcodes.basic[,i]))
      colnames(BetaTable) <- c("Estimate", " ", "Std.Error", " ", "lower", "upper", " ", "p-value", " ")
 
#      significant <- sign(apply(limits,1,prod))
#      starred <- rep("",p)
#      starred[significant==1] <- "*"
#      limits <- as.data.frame(list(limits, starred))
#      colnames(limits)[3] <- " "
  
      llist.basic[[responses[i]]] <- BetaTable
    }
  }

  
  if (confmethod == "BCA")
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, table.bca=llist.bca, method=object$method, conf=object$conf, digits=digits)
  else if (confmethod == "basic")
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, table.basic=llist.basic, method=object$method, conf=object$conf, digits=digits)
  else
    res <- list(responses=responses, covariates=covariates, Betawstd=Betawstd, Sigma=object$Sigma, table.bca=llist.bca, table.basic=llist.basic, method=object$method, conf=object$conf, digits=digits)
}
else
    res <- list(responses=responses, covariates=covariates, Beta=estBeta, Sigma=object$Sigma, method=object$method, digits=digits)

class(res) <- "summary.FRBmultireg"

res

}

###############################################################################

print.summary.FRBmultireg <- function(x, ...) {

if (x$method$est=="MM")
cat(paste("Multivariate regression based on MM-estimates (bdp = ", x$method$bdp, ", eff = ", x$method$eff, ")", sep=""), "\n\n")
else
cat(paste("Multivariate regression based on ", x$method$est, "-estimates (breakdown point = ", x$method$bdp, ")", sep=""), "\n\n")

cat("Response variables: ", x$responses, "\n")
cat("Covariates: ", x$covariates, "\n\n")

if (!is.null(x$Betawstd)) {
  cat("Coefficient estimates (with bootstrap standard errors):\n")
  print(x$Betawstd, digits=x$digits)
} 
else {
  cat("Coefficient estimates):\n")
  print(x$Beta, digits=x$digits)
} 
cat("\nError covariance matrix estimate:\n")
print(x$Sigma, digits=x$digits)
if (!is.null(x$table.bca)) {
    cat("\nConfidence limits (",x$conf*100, "%) and p-values based on BCA method:\n", sep="")
    for (i in 1:length(x$table.bca)) {
    #cat("\n", x$conf*100, "% BCa confidence limits for response \"", x$responses[i], "\":\n\n", sep="")
    cat("\n-response \"", x$responses[i], "\":\n", sep="")
    
    print(x$table.bca[[x$responses[i]]], digits = x$digits)
  }
}

if (!is.null(x$table.basic)) {
    cat("\nConfidence limits (",x$conf*100, "%) and p-values based on Basic Bootstrap:\n", sep="")  
    for (i in 1:length(x$table.basic)) {
    #cat("\n", x$conf*100, "% \"basic bootstrap\" confidence limits for response \"", x$responses[i], "\":\n\n", sep="")
    cat("\n-response \"", x$responses[i], "\":\n", sep="")
    
    print(x$table.basic[[x$responses[i]]], digits = x$digits)
  }
}
cat("\nSignif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")

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