\name{hotelling}

\alias{hotelling}

\title{Simultaneous confidence bounds for a multivariate normal mean using a Hotelling-type confidence ellipse}

\description{xxx}

\usage{hotelling(dat, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{An integer giving the number of steps to use per dimension in the grid search algorithm.}
}

\details{xxx}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{Harold Hotelling (1931) The generalization of Student's ratio. The Annals of Mathematical Statistics, 2(3), 360--378.

Weizhen Wang, J. T. Gene Hwang & Anirban DasGupta (1999) Statistical tests for multivariate bioequivalence. Biometrika, 86(2), 395--402.}

\author{Philip Pallmann}

\note{xxx}

\section{Warning}{Functionality for 3 or more parameters has not yet been thoroughly tested.}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 3), sigma=diag(3) * 0.05)
colnames(mydata) <- c("AUC", "Cmax", "xyz")
hotelling(mydata, steps=10)
}

\keyword{multivariate}
