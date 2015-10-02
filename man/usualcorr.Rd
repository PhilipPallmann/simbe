\name{usualcorr}

\alias{usualcorr}

\title{Simultaneous confidence bounds for a multivariate (?) Gaussian}

\description{xxx}

\usage{usualcorr(dat, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{An integer giving the number of steps to use per dimension in the grid search algorithm.}
}

\details{xxx}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{xxx}

\author{Philip Pallmann}

\note{xxx}

\section{Warning}{Functionality for 3 or more parameters has not yet been throughly tested.}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 3), sigma=(diag(3) * 0.05) + 0.05)
colnames(mydata) <- c("AUC", "Cmax", "xyz")
usualcorr(mydata, steps=10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
