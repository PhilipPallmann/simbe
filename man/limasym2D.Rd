\name{limasym2D}

\alias{limasym2D}

\title{Plot Limacon-shaped bivariate simultaneous confidence region of Brown et al. (1995)}

\description{xxx}

\usage{limasym2D(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3),
          axnames=NULL, main="Limacon (asympt.)")}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
  \item{main}{xxx}
}

\details{xxx}

\value{xxx}

\references{Roger L. Berger & Jason C. Hsu (1996) Bioequivalence trials, intersection-union tests and equivalence confidence sets. Statistical Science, 11(4), 283--319.

Lawrence D. Brown, George Casella and J. T. Gene Hwang (1995) Optimal confidence sets, bioequivalence, and the Limacon of Pascal. Journal of the American Statistical Association, 90(431), 880--889.}

\author{Philip Pallmann}

\note{xxx
%%  ~~further notes~~
}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=diag(2) * 0.05)
limasym2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
          expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
