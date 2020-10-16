#'
#' @description Produce trace and Auto-Correlation Function
#'   plots of MCMC samples
#'  of randomly selected element of A or P or Sigma.
#' @title Trace and ACF plots of randomly selected elements of a variable
#'  in \code{bmrm} object
#' @usage trace_ACF_plot(x,var="P",nplot=3,saveFile=FALSE,...)
#' @param x an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param var name of a variable to which the plots apply.  It sould be one of
#' "A" (source contribution matrix),
#' "P" (source composition or profile matrix),
#' "Sigma" (error variance).
#' @param nplot number of elements to be selected randomly (default=3)
#' @param saveFile TRUE/FALSE, save the plots in file
#' \emph{'var'-trace.pdf} (default=FALSE)
#' @param ... arguments to be passed to methods
#' @details produces trace and auto-correlation function plots of
#' the MCMC samples of randomly selected \code{'nplot'} elements of \code{'var'}
#' in \code{bmrm} object \code{'x'}. If \code{'nplot'}= total number of elements of
#' \code{'var'} then the trace and ACF plots of all elements of 'var' will be provided.
#' @export
#' @examples
#' \dontrun{
#' data(Elpaso)
#' Y=Elpaso$Y ; muP=Elpaso$muP
#' q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#' trace_ACF_plot(out.Elpaso,"Sigma")
#' trace_ACF_plot(out.Elpaso,"P", ACF=T, saveFile=TRUE)
#' trace_ACF_plot(out.Elpaso,"A", nplot=16)
#' }
#'

trace_ACF_plot <- function(x,var="P", ACF=FALSE, nplot=0,saveFile=FALSE,...){

  if (var== "P"  & nplot == 0) nplot=q*ncol( x$Y)
  if (var== "Sigma"  & nplot == 0 ) nplot=ncol( x$Y)
  if (var== "A"  & nplot == 0 ) nplot=16


  var.list<-coda::varnames(x$codaSamples)
  var.list1<-unlist(lapply(var.list,function(x) strsplit(x,"\\[")[[1]][1]))
  id.list<-which(var.list1==var)
  j<-0
  if(!saveFile){
    if(length(id.list)>nplot){
       sel.id<-sample(id.list,size=nplot,replace=FALSE)
     } else{
       sel.id<-id.list
    }
    #grid<-ceiling(sqrt(nplot))

    for(i in sel.id){
      j<-j+1
      if(j %%16 ==1){
        grDevices::X11();
        graphics::par(mfrow=c(4,4))
      }
      y <- coda::mcmc.list(x$codaSamples[,i])
      xp <- as.vector(stats::time(y))
      yp <- if (coda::nvar(y) > 1) {
        y[, j, drop = TRUE]
      } else {
        y
      }
      yp <- do.call("cbind", yp)
      graphics::matplot(xp, yp, xlab = "Iterations", ylab = "", type = 'l',
              col = 4:4,main=var.list[i])
      #ESS<-coda::effectiveSize(y)
      if (ACF==T) stats:: acf(as.matrix(y),main="") #paste0("ESS=",round(ESS,2)))
    }
  } else{
    grDevices::pdf(paste0(var,"-trace_ACF.pdf")) # ,width=6,height=4,paper='special')
    graphics::par(mfrow=c(4,4))

    if(length(id.list)>nplot){
      sel.id<-sample(id.list,size=nplot,replace=FALSE)
    } else{
      sel.id<-id.list
    }
    #grid<-ceiling(sqrt(nplot))

    for(i in sel.id){
      j<-j+1
      y <- coda::mcmc.list(x$codaSamples[,i])
      xp <- as.vector(stats::time(y))
      yp <- if (coda::nvar(y) > 1) {
        y[, j, drop = TRUE]
      } else {
        y
      }
      yp <- do.call("cbind", yp)
      graphics::matplot(xp, yp, xlab = "Iterations", ylab = "", type = 'l',
                        col = 4:4,main=var.list[i])
      #ESS<-coda::effectiveSize(y)
      if(ACF == T) stats:: acf(as.matrix(y),main="") #paste0("ESS=",round(ESS,2)))
       }
    grDevices::dev.off()
    print(paste0("Save as ", getwd(),"/",var,"-trace_ACF.pdf"))
  }
 }
