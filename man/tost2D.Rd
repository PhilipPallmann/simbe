\name{tost2D}

\alias{tost2D}

\title{Plot bivariate (simultaneous) confidence intervals corresponding to two one-sided tests for two parameters}

\description{xxx}

\usage{tost2D(dat, alpha=0.1, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="TOST")}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
  \item{main}{xxx}
}

\details{xxx}

\value{xxx}

\references{Donald Schuirmann (1987) A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and Biopharmaceutics, 15(6), 657--680.}

\author{Philip Pallmann}

\note{xxx}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=diag(2) * 0.05)
tost2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
       expression(paste(Delta, " log(Cmax)"))))
}

\keyword{multivariate}
