#' @description Produce trace and Auto-Correlation Function plots of MCMC samples
#'  of a specific element of A, P, \eqn{\Sigma}.
#' @title Trace and ACF plots of an element of a variable in \code{bmrm} object
#' @usage trace_ACF_plot_indiv(x,var="P",sourceID=1,varID=NULL,obsID=1,...)
#' @param x an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param var name of a variable. It sould be one of
#' "A" (source contribution matrix),
#' "P" (source composition or profile matrix),
#' "Sigma" (error variance).
#' @param sourceID source index (default=1)
#' @param varID variable index (default=1)
#' @param obsID observation index (default=1)
#' @param ... arguments to be passed to methods
#' @details Produce trace and Auto-Correlation Function plots of MCMC samples
#'  of a specific element of A or P or \eqn{\Sigma}. The element is specified by
#'  the indices
#'  \code{sourceID}, \code{varID}, and \code{obsID} of variable \code{var}.
#' @export
#' @examples
#' \dontrun{
#' data(Elpaso)
#' Y=Elpaso$Y ; muP=Elpaso$muP
#' q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#' trace_ACF_plot_indiv(out.Elpaso,"P")
#' trace_ACF_plot_indiv(out.Elpaso,"P",sourceID=2, varID=3 )
#' trace_ACF_plot_indiv(out.Elpaso,"A", obsID=10, sourceID=2 )
#' trace_ACF_plot_indiv(out.Elpaso,"Sigma", varID=3)
#' }
#'
trace_ACF_plot_indiv <- function(x,var="P",sourceID=1,varID=1,obsID=1,...){

  varName = colnames(x$Y)[varID]
  time<-j<-NULL
  #varID.list<-coda::varnames(x$codaSamples)
  parID.list<-coda::varnames(x$codaSamples)
  if(var=="P"){
    id<-paste0("P[",sourceID,",",varID,"]")
    id.name<-paste0("Trace of P: source",sourceID,", ",
                    colnames(x$Y)[varID])
  } else if(var=="A"){
     id<-paste0("A[",obsID,",",sourceID,"]")
     id.name<-paste0("Trace of A: source",sourceID,", observation ID ",obsID)
  } else if(var=="Sigma"){
     id<-paste0("Sigma[",varID,"]")
     id.name<-paste0("Trace of Sigma: ", colnames(x$Y)[varID])
  }

    sel.id<-which(parID.list==id)

  y <- coda::mcmc.list(x$codaSamples[,sel.id])
  xp <- as.vector(stats::time(y))
  yp <- y
  yp <- do.call("cbind", yp)
  graphics::par(mfrow=c(1,2))
  graphics::matplot(xp, yp, xlab = "Iterations", ylab = id, type = 'l',
  #        col = 1:x$nChains,main=id.name)
           col = 4:4,main=id)
  #ESS<-coda::effectiveSize(y)
  stats::acf(as.matrix(y),main=" " ) #paste0("ESS=",round(ESS,2)))

  }

