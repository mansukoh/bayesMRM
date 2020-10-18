#'
#' @description An S3 method that summarizes the
#' output of the \code{bmrm} function in an object of class \code{bmrm}. This object
#'  contains
#' the posterior mean, the posterior standard deviation, and
#' (0.025, 0.05,0.25, 0.5, 0.75,0.95, 0.975) posterior quantiles
#'  of  A, P, \eqn{\Sigma}. It also contains other
#' relevant information about the MCMC procedure such as the burn-in iterations,
#' the number of MCMC chains, etc.
#' @title Summarize the output of the \code{bmrm} function
#' @param object an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param digits integer indicating the humber of signifiant digits
#' @param ... arguments to be passed to methods
#' @export
#' @aliases summary

summary.bmrm <- function(object,digits=3,...){
  if (class(object) != 'bmrm') { stop("incorrect class of 'out.bmrm' ",call.=FALSE)}

  # MCMC
  if (object$errdist=="norm") cat(paste("\n", "E_t ~  normal distribution","\n"))
  if (object$errdist=="t" ) cat(paste("\n","E_t ~ t distribution with df= ",object$df,"\n"))

  cat("Burn-in Iterations = " , object$nBurnIn , ":", object$nBurnIn+object$nIter, "\n")
  cat("Thinning interval = ", object$nThin, "\n")
  cat("Number of MCMC  iterations = ", object$nIter, "\n")

  cat("\n\n")




  # Estimated P
  cat("Estimated Source Composition Profiles 'P':\n\n")
    print(signif(object$P.hat,digits=digits))
  cat("\n\n")

  # Estimated A
  cat("Estimated Source Contribution 'A':\n\n")
  print(utils::head(signif(object$A.hat,digits=digits), n=10))
  cat(" ... with", nrow(object$A.hat)-10, "more rows")
  #as_tibble(rowname_to_column(as.data.frame(object$A.hat)))
  cat("\n\n")

  cat("Estimated Error Variance 'Simga':\n\n")
  print(signif(object$Sigma.hat,digits=digits))
  cat("\n\n")

  cat("Estimated Quantiles of Source Composition Profiles 'P' :\n\n")
  print(utils::head(signif(object$P.quantiles,digits=digits), n=10))
  cat(" ... with", nrow(object$P.quantiles)-10, "more rows")
  cat("\n\n")

  cat("Estimated Quantiles of Source Contribution 'A':\n\n")
  print(utils::head(signif(object$A.quantiles,digits=digits), n=10))
  cat(" ... with", nrow(object$A.quantiles)-10, "more rows")
  cat("\n\n")

  cat("Estimated Quantiles of  Error Variance 'Sigma':\n\n")
  print(signif(object$Sigma.quantiles,digits=digits))
  cat("\n\n")


}

