#' @name bmrm
#' @description  Generate posterior samples of the source
#'  composition matrix  P, the source contribution matrix A,
#'  and the error variance \eqn{\Sigma} using JAGS.
#' @title Bayesian Analysis of Multivariate Receptor Modeling
#' @usage bmrm(Y, q, muP,errdist="norm", df=4,
#'             varP.free=100, xi=NULL, Omega=NULL,
#'               a0=0.01, b0=0.01,
#'              nAdapt=1000, nBurnIn=5000, nIter=5000, nThin=1,
#'              P.init=NULL, A.init=NULL, Sigma.init=NULL,...)
#' @param Y  data matrix
#' @param q  number of sources. It must be a positive integer.
#' @param muP (q,ncol(Y))-dimensional prior mean matrix for the source
#'   composition matrix P, where q is the number of sources. Each row of muP
#'    should have at least q-1 zero elements
#'    (the identifiability condition C1 in Park and Oh (2015)). For the remaining free elements,
#'    nonnegative numbers can be assigned (default=0.5).
#' @param errdist error distribution: either "norm" for normal distribution or "t"
#' for t distribution (default="norm")
#' @param df degrees of freedom of a t-distribution when errdist="t" (default=4)
#' @param varP.free  scalar value of the prior variance of the free (nonzero) elements
#'   of the source composition matrix P (default=100)
#' @param xi prior mean vector of the q-dimensional source contribution
#' vector at time t  (default=vector of 1's)
#' @param Omega diagonal matrix of the prior variance of the q-dimensional
#' source contribution vector at time t (default=identity matrix)
#' @param a0 shape parameter of the Inverse Gamma prior of the error variance
#'    (default=0.01)
#' @param b0 scale parameter of the Inverse Gamma prior of the error variance
#'    (default=0.01)
#' @param nAdapt  number of iterations for adaptation in JAGS (default=1000)
#' @param nBurnIn number of iterations for the burn-in period in MCMC (default=5000)
#' @param nIter number of iterations per chain for monitoring samples from MCMC
#'   (default=5000). \code{nIter} samples are saved in each chain of M CMC.
#' @param nThin thinning interval for monitoring samples from MCMC (default=1)
#' @param P.init initial value of the source composition matrix P. It should have
#'    zero elements in the same zero positions in muP  and it should
#'    satisfy the identifiability conditions C1-C2 of Park and Oh (2015).
#'    If omitted, the nonzero elements of P.init will be randomly generated
#'    from a uniform distrbution.
#' @param A.init initial value of the source contribution matrix A.
#'    If omitted, it will be calculated from Y and P.init.
#' @param Sigma.init initial value of the error variance.
#'    If omitted, it will be calculated from Y, A.init and P.init.
#' @param ... arguments to be passed to methods
#' @return in \code{bmrm} object
#' \describe{
#'   \item{nsource}{number of sources}
#'   \item{nobs}{number of observations in data Y}
#'   \item{nvar}{number of variables in data Y}
#'   \item{Y}{observed data matrix}
#'   \item{errdist}{error distribution}
#'   \item{df}{degrees of freedom when errdist="t"}
#'   \item{A.hat}{posterior mean of the source contribution matrix A}
#'   \item{P.hat}{posterior mean of the source composition matrix P}
#'   \item{Sigma.hat}{posterior mean of the error variance Sigma}
#'   \item{A.sd}{posterior standard deviation of the source contribution matrix A}
#'   \item{P.sd}{posterior standard deviation of the source composition matrix P}
#'   \item{Sigma.sd}{posterior standard deviation of the error variance Sigma}
#'   \item{A.quantiles}{posterior quantlies of A for prob=(0.025, 0.05,
#'   0.25, 0.5, 0.75, 0.95, 0.975)}
#'   \item{P.quantiles}{posterior quantiles of P for prob=(0.025, 0.05,
#'   0.25, 0.5, 0.75, 0.95, 0.975)}
#'   \item{Sigma.quantiles}{posterior quantiles of Sigma for prob=(0.025, 0.05,
#'   0.25, 0.5, 0.75, 0.95, 0.975)}
#'   \item{Y.hat}{estimate of Y computed from A.hat*P.hat}
#'   \item{residual}{Y-Y.hat}
#'   \item{codaSamples}{MCMC posterior samples of A, P, and \eqn{\Sigma} in class "mcmc.list"}
#'   \item{nIter}{number of MCMC iterations per chain for monitoring samples from MCMC}
#'   \item{nBurnIn}{number of iterations for the burn-in period in MCMC}
#'   \item{nThin}{thinning interval for monitoring samples from MCMC}
#' }
#'
#'
#' @details
#' \emph{Model}
#' The basic model for Bayesian multivariate receptor model is
#'   as follows:
#'
#'        \eqn{Y_t=A_t P+E_t, t=1,\cdots,T},
#'
#'   where
#'   \itemize{
#'     \item \eqn{Y_t} is a vector of observations of \eqn{J} variables at time
#'          \eqn{t}, \eqn{t = 1,\cdots,T}.
#'     \item \eqn{P} is a \eqn{q \times J}  source composition
#'             matrix in which the \eqn{k}-th row represents the \eqn{k}-th source
#'            composition profiles, \eqn{k=1,\cdots,q}, \eqn{q} is the number of sources.
#'     \item \eqn{A_t} is a J dimensional source contribution vector at time \eqn{t},
#'            \eqn{t=1,\cdots,T}.
#'     \item \eqn{E_t =(E_{t1}, \cdots, E_{tJ})} is an error term
#'           for the \eqn{t}-th observations,
#'           following \eqn{E_{t} \sim N(0, \Sigma)} or \eqn{E_{t} \sim t_{df}(0, \Sigma)},
#'          independently for \eqn{j = 1,\cdots,J}, where \eqn{\Sigma = diag(\sigma_{1},...,
#'      \sigma_{j})}.
#'  }
#'
#' \emph{Priors}
#'  \itemize{
#'   \item Prior distribution of \eqn{A_t} is given as
#'   a truncated multivariate normal distribution,
#'        \itemize{
#'         \item \eqn{ A_t \sim N(\xi,\Omega) I(A_t \ge 0)}, independently for
#'        \eqn{t = 1,\cdots,T}.
#'        }
#'   \item Prior distribution of \eqn{P_{kj}} (the \eqn{(k,j)}-th element  of the source
#'   composition matrix \eqn{P}) is given as
#'       \itemize{
#'       \item
#'       \eqn{ P_{kj} \sim N(\code{muP}_{kj} , \code{varP.free} )I(P_{kj}
#'       \ge 0)}, for free (nonzero) \eqn{P_{kj}},
#'       \item
#'       \eqn{ P_{kj} \sim N(0, 1e-10 )I(P_{kj} \ge 0)}, for zero \eqn{P_{kj}},
#'
#'       independently for  \eqn{k = 1,\cdots,q; j = 1,\cdots,J }.
#'       }
#'   \item Prior distribution of \eqn{\sigma_j} is \eqn{IG(a0, b0)}, i.e.,
#'      \itemize{
#'       \item \eqn{1/\sigma_j} is \eqn{Gamma(a0, b0)} having mean \eqn{a0/b0},
#'   independently for \eqn{j=1,...,J}.
#'      }
#' }
#'
#' \emph{Notes}
#'   \itemize{
#'       \item
#'      We use the prior
#'      \eqn{ P_{kj} \sim N(0, 1e-10 )I(P_{kj} \ge 0)}  that is practically equal
#'      to the point mass at 0 to simplify the model building in JAGS.
#'      \item  The MCMC samples of A and P are post-processed (rescaled) before saving
#'      so that \eqn{ \sum_{j=1}^J P_{kj} =1} for each \eqn{k=1,...,q} (the identifiablity
#'      condition C3 of Park and Oh (2015).
#'  }
#'
#' @references Park, E.S. and Oh, M-S. (2015), Robust Bayesian Multivariate
#'  Receptor Modeling, CHEMOMETRICS AND INTELLIGENT LABORATORY SYSTEMS,
#'  149, 215-226.
#' @references Park, E.S. and Oh, M-S. (2016), Bayesian Quantile Multivariate
#'  Receptor Modeling, CHEMOMETRICS AND INTELLIGENT LABORATORY SYSTEMS,
#'  159, 174-180.
#'
#' @references Plummer, M. 2003. JAGS: A program for analysis of Bayesian
#' graphical models using Gibbs sampling. Proceedings of the 3rd
#' international workshop on distributed statistical computing, pp. 125.
#'  Technische Universit at Wien, Wien, Austria.
#' @references Plummer, M. 2015. JAGS Version 4.0.0 user manual.

#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data(Elpaso); Y=Elpaso$Y ; muP=Elpaso$muP ; q=nrow(muP)
#' out.Elpaso <- bmrm(Y,q,muP, nAdapt=1000,nBurnIn=5000,nIter=5000,nThin=1)
#' summary(out.Elpaso)
#' plot(out.Elpaso)
#'
#' }

bmrm = function(Y,q,muP,errdist="norm", df=4,
                varP.free=100, xi=NULL, Omega=NULL,
                 a0=0.01, b0=0.01,
                 nAdapt=1000, nBurnIn=5000, nIter=5000, nThin=1,
              P.init=NULL, A.init=NULL, Sigma.init=NULL,...){

  T = nrow(Y)
  J = ncol(Y)


  # =====================================================================
  #  1. Default Setting When input values are omitted
  # =====================================================================

  # convert Y and muP as matrix forms if they are not
  if (!('matrix' %in% class(muP))){ muP <- base::as.matrix(muP) }

  # ==================================================================== #
  # 2. identifiability condition check for muP
  # ==================================================================== #
  idCond_1=1
  for (k in 1:q){
    idCol = which(muP[k,]==0)
    idCond_1 = idCond_1*(length(idCol)  >= q-1)
  }
  if (idCond_1 != 1){ stop("muP violates the identifiability condition C1
  (There are at least   q-1 zero elements in each row of P).") }


  #================================
  if (!('matrix' %in% class(Y))){ Y <- base::as.matrix(Y) }

  # check dimensions of Y and muP

  if ( base::ncol(muP) != J){ stop("incorrect dimension of muP", call.=FALSE) }

  # if xi is missing
  if ( base::is.null(xi) == TRUE ){
    base::cat("Warning: 'xi' is missing. 'xi' is set to default (1).", "\n")
    xi = base::rep(1, q)
  }

  # if Omega is missing
  if ( is.null(Omega) == TRUE ){
    base::cat("Warning: 'Omega' is missing.'Omega' is set to default (1).","\n")
    Omega = base::diag(1,q)
  }

  # if errdist is neither "norm" nor "t"
  if ( errdist != "norm" & errdist != "t" ){
    base::stop("'errdist' is missing or not correct", call.=FALSE)
  }

  # invSigP= inverse of Var of P
  varP.zero=1e-10
  invSigP = base::matrix(1/varP.free, q, J)
  invSigP[which(muP==0)] = 1/varP.zero

  # ======================================================================== #
  #  initial values of P and A
  # ======================================================================== #

  if( base::is.null(P.init) == TRUE ){
    P.init = matrix(stats::runif(q*J),q,J)
    P.init[which(muP==0)] = 0
  }

  #P.init<-t(apply(P.init,1,function(x) x/sum(x)))
  check=idCond_check(P.init)
  if( check != 1) stop("inital value of P does not satisfy the identifiablity
                       conditions ", call.=FALSE)

  if( base::is.null(A.init)== TRUE ){
    A.init = Y %*%t(P.init) %*% solve(P.init%*% t(P.init))
    A.init = A.init*(A.init>0)+0.00001*(A.init<=0)   # to make A.init > 0
  }
  if(is.null(Sigma.init))
     Sigma.init= as.vector(diag(t(Y-A.init%*%P.init)%*%(Y-A.init%*%P.init)/T))


  # ===================================================================== #
  # 3. MCMC sampling with JAGS
  # ===================================================================== #


  if (errdist == "norm"){

    modelString="
    model{
    for (t in 1:T){
    for (j in 1:J){
    Y[t,j] ~ dnorm(mu[t,j], tau[j])
    mu[t,j] <- inprod(A.star[t,1:q], P.star[1:q,j])
    }
    }

    for (j in 1:J){
    tau[j] ~ dgamma(a0,b0)
    Sigma[j]=1/tau[j]
    }

    for (t in 1:T){
    for(k in 1:q){
    A.star[t,k] ~ dnorm(xi[k],1/Omega[k,k])T(0,)
    }
    }

    for (k in 1:q){
    for (j in 1:J){
    P.star[k,j] ~ dnorm(muP[k,j], invSigP[k,j])T(0,)
    }
    }


   #-- post processing ---
   for( k in 1:q){  P.star.sum[k]<- sum(P.star[k,]) }

   for(k in 1:q){
   for( j in 1:J){
      P[k,j]<- P.star[k,j]/P.star.sum[k]
   }}

   for(t in 1:T){
     for(k in 1:q){
      A[t,k]<- A.star[t,k]*P.star.sum[k]
   }}
    }
    "
    dataList = base::list(T=T, J=J, q=q, Y=Y, muP=muP, invSigP=invSigP, xi=xi,
                          Omega=Omega,a0=a0, b0=b0)
    initsList = base::list(A.star=A.init, P.star=P.init, tau=1/Sigma.init)
  }

  else if (errdist == "t"){

    modelString="
    model{

    for (t in 1:T){
    for (j in 1:J){
    Y[t,j] ~ dt(mu[t,j], tau[j], df)
    mu[t,j] <- inprod(A.star[t,1:q], P.star[1:q,j])
    }
    }

    for (j in 1:J){
    tau[j] ~ dgamma(a0,b0)
    Sigma[j]=1/tau[j]
    }

    for (t in 1:T){
    for(k in 1:q){
    A.star[t,k] ~ dnorm(xi[k], 1/Omega[k,k])T(0,)
    }
    }

    for (k in 1:q){
    for (j in 1:J){
    P.star[k,j] ~ dnorm(muP[k,j], invSigP[k,j])T(0,)
    }
    }

#-- post processing ---
   for( k in 1:q){  P.star.sum[k]<- sum(P.star[k,]) }

   for(k in 1:q){
   for( j in 1:J){
      P[k,j]<- P.star[k,j]/P.star.sum[k]
   }}

   for(t in 1:T){
     for(k in 1:q){
      A[t,k]<- A.star[t,k]*P.star.sum[k]
     }}

    }
    "

    dataList = base::list(T=T, J=J, q=q, Y=Y, muP=muP,df=df, invSigP=invSigP, xi=xi,
                          Omega=Omega,  a0=a0, b0=b0)
    initsList = base::list(A.star=A.init, P.star=P.init, tau=1/Sigma.init)

  }

  jagsModel = rjags::jags.model(textConnection(modelString),
                                data=dataList, inits=initsList,
                                n.chains=1, n.adapt=nAdapt)

  stats::update(jagsModel, n.iter=nBurnIn)
  codaSamples = rjags::coda.samples(jagsModel, variable.names=c("A", "P",
                        "Sigma"), n.iter=nIter, thin=nThin)

  # =================================================================== #
  # Posterior estimates of A and P, Sigma
  # =================================================================== #

  mcmcSamples = base::as.matrix(codaSamples)

  A.samples = mcmcSamples[,1:(T*q)]

  A.mean.vec = base::apply(A.samples,2,mean)
  A.sd.vec = base::apply(A.samples,2,stats::sd)
  A.quantiles.vec = base::apply(A.samples,2,
                                stats::quantile,c(0.025,0.05, 0.25, 0.5,
                                                  0.75, 0.95, 0.975))
  A.hat = base::matrix(A.mean.vec, T, q)
  A.sd = base::matrix(A.sd.vec, T, q)
  A.quantiles =t( base::as.matrix(A.quantiles.vec))

  rownames(A.hat)=rownames(A.sd)=paste0("A[",1:T,",]")
  colnames(A.hat)=colnames(A.sd)=paste0("A[,",1:q,"]")
  rownames(A.quantiles)=paste0("A[",rep(1:T,q),",",rep(1:q,each=T),"]")

  colnames(A.quantiles) = c("2.5%","5%", "25%", "50%", "75%", "95%","97.5%")
  ##---------- Phat -------------
  P.samples = mcmcSamples[,(T*q+1):(T*q+q*J)]

  P.mean.vec = base::apply(P.samples,2,mean)
  P.sd.vec = base::apply(P.samples,2,stats::sd)
  P.quantiles.vec = base::apply(P.samples,2,
                                stats::quantile,c(0.025,0.05, 0.25, 0.5,
                                                  0.75, 0.95, 0.975))
  P.hat = base::matrix(P.mean.vec, q, J)
  P.sd = base::matrix(P.sd.vec, q, J)
  P.quantiles = t(base::as.matrix(P.quantiles.vec))

  rownames(P.hat)=rownames(P.sd)=paste0("P[",1:q,",]")
  colnames(P.hat)=colnames(P.sd)=paste0("P[,",1:J,"]")
  rownames(P.quantiles)=paste0("P[",rep(1:q,J),",",rep(1:J,each=q),"]")
  colnames(P.quantiles) = c("2.5%","5%", "25%", "50%", "75%",
                            "95%","97.5%")


  # Rescale P.hat, P.sd, A.hat, A.sd so that row sum of P.hat == 1
  #rowsum.Phat=apply(P.hat,1,sum)
  #for( k in 1:q){
  #  P.hat[k,]=P.hat[k,]/rowsum.Phat[k]
  #  P.sd[k,]=P.sd[k,]/rowsum.Phat[k]
  #  A.hat[,k]=A.hat[,k]*rowsum.Phat[k]
  #  A.sd[,k]=A.sd[,k]*rowsum.Phat[k]
  #}


  P.hat[which(muP==0)]=0
  check=idCond_check(P.hat)
  if(check != 1) cat("Warning: P.hat does not satisfy the identifiability
                     conditions")

  ##---------- Sigma hat -------------
  Sigma.samples = mcmcSamples[,(T*q+q*J+1):(T*q+q*J+J)]
  Sigma.hat = base::apply(Sigma.samples,2,mean)
  Sigma.sd = base::apply(Sigma.samples,2,stats::sd)
  Sigma.quantiles.vec = base::apply(Sigma.samples,2,
                                    stats::quantile,c(0.025,0.05, 0.25, 0.5,
                                                      0.75, 0.95, 0.975))
  Sigma.quantiles =t( base::as.matrix(Sigma.quantiles.vec))

  names(Sigma.hat) = paste0("Sigma[",1:J,"]")
  names(Sigma.sd) = paste0("Sigma[",1:J,"]")
  rownames(Sigma.quantiles) = paste0("Sigma[",1:J,"]")
  colnames(Sigma.quantiles) = c("2.5%","5%", "25%", "50%", "75%",
                                "95%","97.5%")

  Y.hat = A.hat %*% P.hat
  residual=Y-Y.hat


  # output list
  out<-list()
  out$nsource=q
  out$nobs= T
  out$nvar= J
  out$Y = Y
  out$errdist = errdist
  out$df=df
  out$A.hat = A.hat
  out$P.hat = P.hat
  out$Sigma.hat = Sigma.hat
  out$A.sd = A.sd
  out$P.sd = P.sd
  out$Sigma.sd = Sigma.sd
  out$A.quantiles =A.quantiles
  out$P.quantiles =P.quantiles
  out$Sigma.quantiles =Sigma.quantiles
  out$Y.hat = Y.hat
  out$residual=residual
  out$codaSamples =codaSamples

  out$nIter = nIter
  out$nBurnIn=nBurnIn
  out$nThin=nThin

  base::class(out) = "bmrm"
  return(out)
}


