\name{Elpaso}
\alias{Elpaso}
\docType{data}
\title{
PM2.5 speciation data from El Paso, Texas, USA.
}
\description{
The data frame ahs the following components:

\itemize{
\item \emph{Y} 224 by 15 matrix of 224 observations on 15 PM2.5 species.
PM 2.5 was measured every three days during the time period of
1/2/2006 ~ 4/7/2009 from the Chamizal station in the city of El Paso, USA.
Out of the 58 original PM 2.5 species, 15 species were selected. After removing
any observations with missig values, the final data consists of 224 complete
observatins on the following 15 PM2.5 species.
  \describe{
    \item{Al}{Aluminum}
    \item{Ca}{Calcium}
    \item{Cl2}{Chlorine}
    \item{EC}{EC CSN}
    \item{Fe}{Iron}
    \item{K_p}{Potassium ion}
    \item{Mg}{Magnanese}
    \item{NV_NO3}{Non-volatile nitrate}
    \item{NH4_p}{Ammonium ion}
    \item{Na}{Sodium}
    \item{OC}{OC CSN unadjusted}
    \item{SO4}{Sulfate}
    \item{Si}{Silicon}
    \item{Ti}{Titanium}
    \item{Zn}{Zinc}
}

\item \emph{muP}  4 by 15 matrix of the prior mean of the source composition matrix P for data.  Zero values  are assigned for some elements of muP to satisfy the
identifiability conditions C1-C2 in Park and Oh (2015). The remaining nonzero elements of muP have value 0.5. Note that the number of sources (the number of rows in muP) is presumed to be 4 here.
}
}

\references{
 Park, E.S. and Oh, M-S. (2016), Bayesian Quantile Multivariate
  Receptor Modeling, CHEMOMETRICS AND INTELLIGENT LABORATORY SYSTEMS,
  159, 174-180.
}
\examples{
data(Elpaso)
Y=Elpaso.Y
muP=Elpaso.muP
}
\keyword{datasets}
