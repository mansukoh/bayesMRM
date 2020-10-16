#' @description Draw principal component plots of data (Y) and
#' source profiles (rows) of the estimated source composition matrix P.hat
#' (and P0 if there is another source composition matrix P0 to compare, e.g.,
#' P0 could be the true P in simulation or P0 could be another estimate of P)
#' @title Principal component plot
#' @usage pcplot(x, P0, G3D=FALSE,...)
#' @param x an object of class \code{bmrm}, the output from the function \code{bmrm}
#' @param P0 estimated value of P (in simulation it can be the true value of P)
#' @param G3D TRUE/FALSE, dynamic 3D plot (default=FALSE)
#' @param ... arguments to be passed to methods
#' @export
#' @examples
#' \dontrun{
#' data(Elpaso)
#' Y=Elpaso$Y ; muP=Elpaso$muP
#' q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#'
#' pcplot(out.Elpaso)
#' pcplot(out.Elpaso,G3D=TRUE)
#' }


pcplot <- function(x,P0=NULL,G3D=FALSE,...){
   Y <- x$Y
   Yn <- t(apply(Y,1,function(x) x/sum(x)) )

   Phat <- x$P.hat
   Pn <- t(apply(Phat,1,function(x) x/sum(x)) )
   Y.svd <- svd(stats::cor(Y))
   Z <- Yn %*%Y.svd$v
   S <- Pn %*%Y.svd$v

   if(G3D){
      G3D.data<-rbind(Z[,1:3],S[,1:3])
      G3D.color<-c(rep("lightblue",nrow(Z)),rep("red",3))
      G3D.pch<-c(rep(16,nrow(Z)),c(2,3,4))
      G3D.text<-paste0("S",1:nrow(S))
      if(!is.null(P0)){
         P0n <- t(apply(P0,1,function(x) x/sum(x)))
         T <- P0n %*%Y.svd$v
         G3D.data<-rbind(G3D.data,T[,1:3])
         G3D.color<-c(G3D.color,rep("blue",3))
         G3D.pch<-c(G3D.pch,c(2,3,4))
         G3D.text<-c(G3D.text,paste0("P",1:nrow(T)))
      }

      rgl::plot3d(G3D.data[,1:3],col=G3D.color,
         xlab="z1",ylab="z2",zlab="z3",radius=0.005,type="s",family=2)
      rgl::text3d(G3D.data[-(1:nrow(Y)),1:3],text=G3D.text,pos=1,font=2)
   } else{
     ggplot.data<-data.frame(Z)
     colnames(ggplot.data)<-paste0("z",1:ncol(ggplot.data))
     z1<-ggplot.data$z1;z2<-ggplot.data$z2;z3<-ggplot.data$z3
     ggplot.label1<-data.frame(S)
     colnames(ggplot.label1)<-paste0("z",1:ncol(ggplot.label1))
     ggplot.label1$label<-paste0("S",1:nrow(ggplot.label1))
     label<-ggplot.label1$label

     P1<-ggplot2::ggplot(ggplot.data,ggplot2::aes(z1,z2))+
           ggplot2::geom_point(color="lightblue")+
           ggplot2::geom_point(data=ggplot.label1,
              ggplot2::aes(z1,z2),shape=2,size=2,col="blue")+
           ggplot2::geom_text(data=ggplot.label1,ggplot2::aes(label=label))

     P2<-ggplot2::ggplot(ggplot.data,ggplot2::aes(z1,z3))+
           ggplot2::geom_point(color="lightblue")+
           ggplot2::geom_point(data=ggplot.label1,
              ggplot2::aes(z1,z3),shape=2,size=2,col="blue")+
           ggplot2::geom_text(data=ggplot.label1,ggplot2::aes(label=label))

     P3<-ggplot2::ggplot(ggplot.data,ggplot2::aes(z2,z3))+
           ggplot2::geom_point(color="lightblue")+
           ggplot2::geom_point(data=ggplot.label1,
              ggplot2::aes(z2,z3),shape=2,size=2,col="blue")+
           ggplot2::geom_text(data=ggplot.label1,ggplot2::aes(label=label))

     if(!is.null(P0)){
       P0n <- t(apply(P0,1,function(x) x/sum(x)))
       T <- P0n %*%Y.svd$v
       ggplot.label2<-data.frame(T)
       colnames(ggplot.label2)<-paste0("z",1:ncol(ggplot.label2))
       ggplot.label2$label<-paste0("P",1:nrow(ggplot.label2))
       P1<-P1+ggplot2::geom_text(data=ggplot.label2,ggplot2::aes(label=label))
       P2<-P2+ggplot2::geom_text(data=ggplot.label2,ggplot2::aes(label=label))
       P3<-P3+ggplot2::geom_text(data=ggplot.label2,ggplot2::aes(label=label))
     }
     gridExtra::grid.arrange(P1,P2,P3,ncol=2)
  }
}
