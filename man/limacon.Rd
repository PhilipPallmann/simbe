\name{limacon}

\alias{limacon}

\title{Simultaneous confidence bounds for a multivariate (?) Gaussian}

\description{xxx}

\usage{limacon(dat, alpha=0.1, steps=100)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{An integer giving the number of steps to use per dimension in the grid search algorithm.}
}

\details{xxx}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{Roger L. Berger & Jason C. Hsu (1996) Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science, 11(4), 283--319.

Lawrence D. Brown, George Casella and J. T. Gene Hwang (1995) Optimal confidence sets, bioequivalence, and the Limacon of Pascal. Journal of the American Statistical Association, 90(431), 880--889.}

\author{Philip Pallmann}

\note{xxx}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 3), sigma=diag(3) * 0.05)
colnames(mydata) <- c("AUC", "Cmax", "xyz")
limacon(mydata, steps=10)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
