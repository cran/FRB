
summary.FRBhot <- function(object, digits=5, ...) {
  res <- list(object=object, digits=digits)
  class(res) <- "summary.FRBhot"
  res
}

print.summary.FRBhot <- function(x, ...) {

  digits <- x$digits
  x <- x$object

  if (is.null(x$Mu1))
    {cat(x$meth, "\n\n")
     cat("data: ", deparse(x$data),"\n")
     cat("T^2_R = ",round(x$teststat,2), "\n")
     cat("p-value = ", round(x$pvalue,4), "\n")
#     cat("Alternative hypothesis : true mean vector is not equal to (",round(x$Mu0,6), ")", "\n")
     cat("Alternative hypothesis : true mean vector is not equal to (", x$Mu0, ")", "\n")
     cat("\n", x$conf*100, "% simultaneous confidence intervals for components of mean :\n")
     print(x$CI,digits=digits)
     cat("\nSample estimates :\n")
     cat("   location:\n")
     print(x$Mu, digits=digits)
     cat("\n   covariance:\n")
     print(x$Sigma, digits=3)
  }
  else { cat(x$meth, "\n\n")
       cat("data: ", deparse(x$data[[1]]), " and ", deparse(x$data[[2]]),"\n")
       cat("T^2_R = ",round(x$teststat,2),  "\n")
       cat("p-value = ", round(x$pvalue,4), "\n")
       cat("Alternative hypothesis : true difference in mean vectors is not equal to (",rep(0,length(x$Mu1)), ")", "\n")
       cat("\n", x$conf*100, "% simultaneous confidence intervals for components of difference in means :\n")
       print(x$CI,digits=3)
       cat("\nSample estimates :\n")
       cat("   location:\n")
       printmat <- rbind(x$Mu1,x$Mu2,x$Mu1-x$Mu2)
       rownames(printmat) <- c("   Sample 1","   Sample 2", " difference")
       print(printmat, digits=digits)
       cat("\n   covariance:\n")
       print(x$Sigma, digits=3)
    }
}

print.FRBhot <- function(x, digits=5, ...)  {

if (is.null(x$Mu1))
    {cat(x$meth, "\n\n")
     cat("data: ", deparse(x$data),"\n")
     cat("T^2_R = ",round(x$teststat,2), "\n")
     cat("p-value = ", round(x$pvalue,4), "\n")
     cat("Alternative hypothesis : true mean vector is not equal to (", x$Mu0, ")", "\n")
}
else { cat(x$meth, "\n\n")
       cat("data: ", deparse(x$data[[1]]), " and ", deparse(x$data[[2]]),"\n")
       cat("T^2_R = ",round(x$teststat,2),  "\n")
       cat("p-value = ", round(x$pvalue,4), "\n")
       cat("Alternative hypothesis : true difference in mean vectors is not equal to (",rep(0,length(x$Mu1)), ")", "\n")
}

}
