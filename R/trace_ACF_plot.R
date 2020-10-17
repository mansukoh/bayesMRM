#'
#' @description Produce trace and Auto-Correlation Function
#'   plots of MCMC samples
#'  of  elements of A or P or Sigma.
#' @title Trace and/or ACF plots of elements of a variable
#'  in \code{bmrm} object
#' @usage trace_ACF_plot(x,var="P",ACF=FALSE, nplot=16,saveFile=FALSE,...)
#' @param x an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param var name of a variable to which the plots apply.  It sould be one of
#' "A" (source contribution matrix),
#' "P" (source composition or profile matrix),
#' "Sigma" (error variance).
#' @param ACF TRUE/FALSE  IF TRUE ACF plot will be provided along with trace
#'  plot (dafault: FALSE)
#' @param nplot  number of elements of 'var' for trace/ACF plot. If 'nplot' is
#'  smaller than the total number of elements of 'var' then trace/ACF plots of
#'  'nplot' randomly selected elements will be drawn.
#'   Otherwise, trace/ACF plots of  all elements will be drawn.
#'    (default=0  implies that all elements will be selected if var="P" or "Sigma"
#'    and 16 elements will be selected  randomly if  var="A")
#' @param saveFile TRUE/FALSE, save the plots in file
#' \emph{'var'-trace.pdf} (default=FALSE)
#' @param ... arguments to be passed to methods
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

  if (var== "P"  & nplot == 0) nplot=nrow(x$P.hat)*ncol( x$P.hat)
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

    #par("mar")
    par(mar=rep(1,4))


    for(i in sel.id){
      j<-j+1
      if(j %%(4*4) ==1){
         if( j >1 ) grDevices::X11();
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
      graphics::matplot(xp, yp, xlab = "Iteration", ylab = "", type = 'l',
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
      graphics::matplot(xp, yp, xlab = "Iteration", ylab = "", type = 'l',
                        col = 4:4,main=var.list[i])
      #ESS<-coda::effectiveSize(y)
      if(ACF == T) stats:: acf(as.matrix(y),main="") #paste0("ESS=",round(ESS,2)))
       }
    grDevices::dev.off()
    print(paste0("Save as ", getwd(),"/",var,"-trace_ACF.pdf"))
  }
 }
