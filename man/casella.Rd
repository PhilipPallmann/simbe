\name{casella}

\alias{casella}

\title{Simultaneous confidence bounds for a multivariate normal mean using Casella & Hwang's empirical Bayes method}

\description{xxx}

\usage{casella(dat, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{An integer giving the number of steps to use per dimension in the grid search algorithm.}
}

\details{xxx}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{George Casella & Jiunn Tzon Hwang (1983) Empirical Bayes confidence sets for the mean of a multivariate normal distribution. Journal of the American Statistical Association, 78(383), 688--698.}

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
casella(mydata, steps=10)
}

\keyword{multivariate}
