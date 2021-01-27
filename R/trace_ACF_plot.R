#'
#' @description Produce trace and Auto-Correlation Function
#'   plots (with Effective sample size) of MCMC samples
#'  of  elements of A, nonzero elements of P, elements of Sigma.
#' @title Trace and/or ACF plots of elements of a variable
#'  in \code{bmrm} object
#' @usage trace_ACF_plot(x,var="P",ACF=FALSE, nplot=12,irow=2, icol=3,saveFile=FALSE,...)
#' @param x an object of class \code{bmrm}, the output of the \code{bmrm} function
#' @param var name of a variable to which the plots apply.  It should be one of
#' "A" (source contribution matrix),
#' "P" (source composition matrix),
#' "Sigma" (error variance).
#' @param ACF TRUE/FALSE  IF TRUE ACF plot will be provided and Effective sample size(dafault: FALSE)
#' @param nplot  number of elements of 'var' for trace and/or ACF plot. If 'nplot' is
#'  smaller than the total number of elements of 'var' then plots of
#'  'nplot' selected elements will be drawn. Otherwise, trace/ACF plots of  all
#'   elements will be drawn.
#'    (default=0  implies that all elements will be selected if var="P" or "Sigma"
#'    and the first 12 elements will be selected  if  var="A")
#' @param irow row index of A/P matrix or index of element of Sigma vector.
#'     Plots of 'nplot' elements starting from (irow, icol) element of A/P or
#'     elements starting from irow element of Sigma will be drawn (default=1).
#' @param icol column number of  A/P matrix.  Plots of 'nplot' elements starting
#' from (irow, icol) element of A/P will be drawn (default=1).
#' @param saveFile TRUE/FALSE, save the plots in file
#' \emph{'var'-trace.pdf} (default=FALSE)
#' @param ... arguments to be passed to methods
#' @export
#' @examples
#' \dontrun{
#' data(Elpaso); Y=Elpaso$Y ; muP=Elpaso$muP ; q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#' trace_ACF_plot(out.Elpaso,"Sigma", ACF=T)
#' trace_ACF_plot(out.Elpaso,"P", ACF=T)
#' trace_ACF_plot(out.Elpaso,"P", ACF=T,saveFile=T )
#' trace_ACF_plot(out.Elpaso,"A",ACF=T, nplot=12, irow=2, icol=3)
#' }
#'

trace_ACF_plot <- function(x,var="P", ACF=FALSE, nplot=0,irow=1, icol=1, saveFile=FALSE,...){

  if (var== "P"  & nplot == 0) nplot=nrow(x$P.hat)*ncol( x$P.hat)
  if (var== "Sigma"  & nplot == 0 ) nplot=ncol( x$Y)
  if (var== "A"  & nplot == 0 ) nplot=12

  var.list<-coda::varnames(x$codaSamples)
  var.list1<-unlist(lapply(var.list,function(x) strsplit(x,"\\[")[[1]][1]))
  id.list<-which(var.list1==var)
  if (var =="P" ){ id.list= id.list[as.vector(x$muP) != 0  ]   }

  nplot=min( nplot, length(id.list))

  istart=1
  if( var =="P" & nplot>0 & nplot< length(id.list)) {
       istart= (icol-1)*x$nsource + irow
  }
  if( var =="A" & nplot>0 & nplot< length(id.list)) {
    istart= (icol-1)*x$nobs + irow
  }
  if( var =="Sigma" & nplot>0 & nplot< length(id.list)) {
    istart= irow
  }

  #if( ACF==T ) nplot=as.integer(nplot/2)

  j<-0
  if(!saveFile){
    if(length(id.list)>nplot){
           sel.id<- id.list[istart:(istart-1+nplot)] #sample(id.list,size=nplot,replace=FALSE)
     } else{
       sel.id<-id.list
    }

    graphics::par(mfrow=c(3,4))


    for(i in sel.id){
      j<-j+1

      if( ACF==F & j %%(3*4) ==1  &  j >1 ){
       grDevices::X11()
       graphics::par(mfrow=c(3,4))
       graphics::par(mar=rep(2,4))
      }

      if( ACF==T & j %%(3*2) ==1 & j >1 ){
        grDevices::X11()
        graphics::par(mfrow=c(3,4))
        graphics::par(mar=rep(2,4))
      }

      y <- coda::mcmc.list(x$codaSamples[,i])
      xp <- as.vector(stats::time(y))
      yp <- if (coda::nvar(y) > 1) {
        y[, j, drop = TRUE]
      } else {
        y
      }
      yp <- do.call("cbind", yp)
      graphics::matplot(xp,yp,xlab="Iteration",ylab="",type ='l',
           col = 4:4,main=var.list[i])
         ESS<-coda::effectiveSize(y)
      if (ACF==T) stats::acf(as.matrix(y),main=paste0("ESS=",round(ESS,2)))
    }

  } else {

    grDevices::pdf(paste0(var,"-trace_ACF.pdf")) #,width=6,height=4,paper='special')
    graphics::par(mfrow=c(3,4))
    graphics::par(mar=rep(2,4))

    if(length(id.list)>nplot){
      sel.id<- id.list[istart:(istart-1+nplot)]
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
      ESS<-coda::effectiveSize(y)
      if(ACF == T) stats:: acf(as.matrix(y),main=paste0("ESS=",round(ESS,2)))
    }

    grDevices::dev.off()
    print(paste0("Save as ", getwd(),"/",var,"-trace_ACF.pdf"))
  }

  }
