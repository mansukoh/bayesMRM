#' Convergence Diagnostics on MCMC samples in \code{bmrm}
#' @description Compute convergence diagnostic measures of 
#'  Geweke (1992), Heidelberger and Welch (1983).
#' @usage convdiag_bmrm(x , var="P", convdiag="heidel",print=TRUE,...)
#' @param x an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param var name of a variable for convergence disagnostics. It should be one of "A" (source contribution matrix), 
#' "P" (source composition or profile matrix), or
#' "Sigma" (error variance).
#' @param convdiag  vector of convergence diagnostic measures. It should be any subvector
#'  of ("geweke", "heidel") (default="heidel").
#' @param print TRUE/FALSE, print convergence diagnostics results (default=TRUE)
#' @param ... arguments to be passed to methods
#' 
#' @return A list of  convergence diagnostics results
#' \describe{
#'   \item{convdiag}{selected convergence diagnostic measures} 
#'   \item{geweke}{Geweke's statistics and p-value if \code{convdiag}
#'  includes "geweke", NULL if \code{convdiag} does not include "geweke"} 
#'   \item{heidel}{Heidelberger and Welch's test statistics 
#' and p-value if \code{convdiag} includes "heidel"; NULL if 
#' \code{convdiag} does not include "heidel"} 
#' }
#'
#' @details     The functions  \code{geweke.diag}, \code{heidel.diag} 
#'  in \bold{coda} package are used to compute these diagnostics. 
#'  In our package, \code{frac1=0.1, frac2=0.5} are used as parameters of the function
#'  \code{coda::geweke.diag} and \code{eps = 0.1, pvalue = 0.05} are used as
#'   parameters of the function \code{coda::heidel.diag}.

#' @references Geweke, J.(1992) Evaluating the accuracy of sampling-based
#' approaches to calculating posterior moments. In Bayesian Statistics 4
#' (ed JM Bernado, JO Berger, AP Dawid and AFM Smith). Clarendon Press.
#' @references Heidelberger P, and Welch PD. (1981) A spectral method for
#' confidence interval generation and run length control in simulations.
#' Comm. ACM. 24, 233-245.
#' @references Heidelberger P. and Welch PD.(1983) Simulation run length
#' control in the presence of an initial transient.
#' Opns Res., 31, 1109-44,Oxford, UK.
#' @references Plummer, M., Best, N., Cowles, K. and Vines K. (2006) CODA: 
#' Convergence Diagnosis and Output Analysis for MCMC, R News, Vol 6, pp. 7-11.
#'
#' @export
#' @examples
#' \dontrun{
#' data(Elpaso)
#' Y=Elpaso$Y ; muP=Elpaso$muP
#' q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#' conv1<-convdiag_bmrm(out.Elpaso,var="P")
#' conv2<-convdiag_bmrm(out.Elpaso,var="A", convdiag="geweke")
#' conv3<-convdiag_bmrm(out.Elpaso,var="Sigma", convdiag=c("geweke","heidel"))
#' }


convdiag_bmrm <- function(x , var="P", convdiag="heidel",print=TRUE,...){
  
  if (class(x) != 'bmrm') { stop("incorrect class of 'x.bmrm' ",call.=FALSE)}
  T<-x$nobs
  q<-x$nsource
  J<-x$nvar
  var.codaSamples=list()
  geweke_table = heidel_table = NULL
  if (var == "P"){
        for (ic in 1:length(x$codaSamples)) {
           var.codaSamples[[ic]]=x$codaSamples[[ic]][,(T*q+1):(T*q+q*J)] }
  } else if (var == "A"){
    for (ic in 1:length(x$codaSamples)) {
      var.codaSamples[[ic]]=x$codaSamples[[ic]][,1:(T*q)] }
  } else if (var == "Sigma"){
      for (ic in 1:length(x$codaSamples)) {
        var.codaSamples[[ic]]=x$codaSamples[[ic]][,(T*q+q*J+1):(T*q+q*J+J)] }
  } else{
      cat("unknown 'var' name \n" )}

  class(var.codaSamples)<-"mcmc.list"
  if (!("geweke" %in% convdiag)& 
      !("heidel" %in% convdiag)) { 
    stop("incorrect type of convergence diagnostics", call.=FALSE) }


  if ("geweke" %in% convdiag){
  #geweke: a test for equality of the means of the first and last part of
  #a Markov chain (by default the first 10% and the last 50%)

    mcmcSamples = base::as.matrix(var.codaSamples)
    geweke1 = coda::geweke.diag(mcmcSamples, frac1=0.1, frac2=0.5)
    geweke_pvalue <- 2*stats::pnorm(-base::abs(geweke1$z))
    geweke_table <- base::data.frame('Geweke statistics' = geweke1$z,
                                     'Geweke p value' =
                                      round(geweke_pvalue, digits = 4))
    if(print){
      cat("\n\n")
      cat("Geweke Diagnostics :\n\n")
      base::print(utils::head(geweke_table,n=10))
      if(nrow(geweke_table)>10)
        base::cat(paste("...",nrow(geweke_table)-10," more rows"))    
      cat("\n\n")
    }
  } 
  
  if ("heidel" %in% convdiag){
    mcmcSamples=base::as.matrix(var.codaSamples)
    
    heidel1 <- coda::heidel.diag(mcmcSamples, eps=0.1, pvalue=0.05)
    heidel_table <- base::data.frame('Heidel test' = heidel1[,1],
                                        'Heidel p value' =round(heidel1[,3], 4))
    if(print){
      base::cat("\n\n")
      base::cat("Heidel Diagnostics:\n\n")
      base::print(utils::head(heidel_table,n=10))
      if(nrow(heidel_table)>10)
        base::cat(paste("...",nrow(heidel_table)-10," more rows"))          
      base::cat("\n\n" )
    }
  }
  
  return(list(convdiag = convdiag,
              geweke=geweke_table,
              heidel=heidel_table))
}

