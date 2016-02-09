\name{standardcorr2D}

\alias{standardcorr2D}

\title{Plot bivariate simultaneous confidence region xxx}

\description{xxx}

\usage{standardcorr2D(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3),
                      axnames=NULL, main="Standard (corr.)")}

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

\references{xxx}

\author{Philip Pallmann}

\note{xxx
%%  ~~further notes~~
}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=(diag(2) * 0.05) + 0.05)
standardcorr2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
               expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
