#'
#' @description Produce plots of the estimated
#' posterior mean and 95\% posterior intervals of A,P, Sigma
#' based on the MCMC samples in \code{bmrm}.
#' @title Produce plots of the parameter estimates
#' @param x an object of class \code{bmrm}, the output of the function \code{bmrm}
#' @param type name of a variable (default="P").
#' It should be one of "P"(source composition or profile matrix P),
#' "A"(source contribution matrix A), "both" (both P and A), "Sigma" (error  variance).
#' @param text TRUE/FALSE, display the value of P.hat on the plot for P (defaut=FALSE)
#' @param ... arguments to be passed to methods
#'
#' @details
#' The following types of plots are drawn depending on the selected parameters:
#' \itemize{
#' \item P: bar plots of the posterior medians with 95\% posterior intervals of elements for each row of P
#' \item A: time series plots of  posterior medians with 95\% posterior intervals elements for each column of A
#' \item Sigma: bar plots of  posterior medians with 95\% posterior intervals of elements of Sigma
#'  }
#' @export
#' @aliases plot
#'


plot.bmrm <- function(x,type="both",...){


  BB<-A<-K<-Pname<-P<-LB<-UB<-Sigma<-NULL
## A
  if(type=="A" | type=="both"){
    ggplot.data<-data.frame(K=paste("source",rep(1:x$nsource,each=nrow(x$Y))),
      T=rep(1:x$nobs,x$nsource),
      LB=x$A.quantiles[,"2.5%"],
      Med=x$A.quantiles[,"50%"],
      UB=x$A.quantiles[,"97.5%"],
      A=c(x$A.hat))

    ggplot.polyg.data<-data.frame(T=c(ggplot.data$T,
                                    ggplot.data$T[nrow(ggplot.data):1]),
                BB=c(ggplot.data$LB,ggplot.data$UB[nrow(ggplot.data):1]),
                K=c(as.character(ggplot.data$K),
                    as.character(ggplot.data$K[nrow(ggplot.data):1])))

    P1<-ggplot2::ggplot(ggplot.data)+
            ggplot2::facet_grid(K~.,scales="free_y")+
            ggplot2::geom_polygon(data=ggplot.polyg.data,
                ggplot2::aes(T,BB),alpha=0.5,col="pink",fill="pink")+
            ggplot2::geom_line(ggplot2::aes(T,A))+
            ggplot2::ylab("A") +
            ggplot2::xlab("obs")
  }
  if(type=="P"| type=="both"){
    ggplot.data<-data.frame(K=paste("source",rep(1:x$nsource,each=x$nvar)),
                            #LB=x$P.quantiles[,"2.5%"]*100,
                            LB=c( t( matrix(x$P.quantiles[,"2.5%"]*100,
                                           nrow(x$P.hat), ncol(x$P.hat))) ) ,
                            #Med=x$P.quantiles[,"50%"]*100,
                            #UB=x$P.quantiles[,"97.5%"]*100,
                            UB=c( t( matrix(x$P.quantiles[,"97.5%"]*100,
                                            nrow(x$P.hat), ncol(x$P.hat))) ) ,
                            P=c(t(x$P.hat)*100),
                            Pname=rep(colnames(x$Y),x$nsource),
                            Pid = rep(1:ncol(x$Y),x$nsource))  #each=x$nsource))

    P2<-ggplot2::ggplot(ggplot.data,ggplot2::aes(factor(Pid),P))+
          ggplot2::geom_bar(stat="identity",width=0.5,fill="gray60")+
          ggplot2::facet_grid(K~.,scales="free_y")+ggplot2::xlab("variables")+
          ggplot2::geom_errorbar(ggplot2::aes(ymin=LB,ymax=UB),
                       color="gray40",width=0.2)+
         ggplot2::scale_x_discrete(labels=colnames(x$Y))
    #if(text)
    #  P2<-P2+ggplot2::geom_text(ggplot2::aes(label=as.character(round(P,2))),
    #              nudge_y=3,nudge_x=-0.3)
  }
  if(type=="Sigma"){
    ggplot.data<-data.frame(LB=x$Sigma.quantiles[,"2.5%"]*100,
                            Med=x$Sigma.quantiles[,"50%"]*100,
                            UB=x$Sigma.quantiles[,"97.5%"]*100,
                            Sigma=c(x$Sigma.hat*100),
                            Pname=rep(colnames(x$Y)))

    P3<-ggplot2::ggplot(ggplot.data,ggplot2::aes(Pname,Sigma))+
          ggplot2::geom_point()+
          ggplot2::xlab("variables")+
          ggplot2::geom_errorbar(ggplot2::aes(ymin=LB,ymax=UB),
                       color="gray40",width=0.2)
  }

  P<-list()
  if(type=="A"){
    P[[1]]<-P1  #print(P1)
  } else if(type=="P"){
    P[[1]]<-P2   #print(P2)
  } else if(type=="Sigma"){
    P[[1]]<-P3   #print(P3)
  } else if(type=="both"){
    P[[1]]<-P1; P[[2]]<-P2  #print(P1);print(P2)
  }
  gridExtra::marrangeGrob(grobs=P,nrow=1,ncol=1,top=NULL)
}

