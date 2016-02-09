\name{tseng}

\alias{tseng}

\title{Simultaneous confidence bounds for a multivariate normal mean using Tseng's minimum expected length method}

\description{xxx}

\usage{tseng(dat, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{An integer giving the number of steps to use per dimension in the grid search algorithm.}
}

\details{xxx}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{Yu-Ling Tseng (2002) Optimal confidence sets for testing average bioequivalence. Test, 11(1), 127--141.}

\author{Philip Pallmann}

\section{Warning}{Functionality for 3 or more parameters has not yet been thoroughly tested.}

\note{When the estimates are close to zero, the confidence region may not exist.}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 3), sigma=diag(3) * 0.05)
colnames(mydata) <- c("AUC", "Cmax", "xyz")
tseng(mydata, steps=10)
}

\keyword{multivariate}
