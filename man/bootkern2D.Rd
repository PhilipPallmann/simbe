\name{bootkern2D}

\alias{bootkern2D}

\title{Plot bivariate simultaneous confidence region using nonparametric bootstrapping and bivariate kernel density estimation}

\description{xxx}

\usage{bootkern2D(dat, alpha=0.1, nboot=10000, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Bootstrap")}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{nboot}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
  \item{main}{xxx}
}

\details{xxx}

\value{xxx}

\references{A. C. Davison and D. V. Hinkley (1997) Bootstrap Methods and their Application. Cambridge University Press, New York, NY. ISBN 978-0-521-57391-2.

Bradley Efron and Robert J. Tibshirani (1993) An Introduction to the Bootstrap. Chapman and Hall/CRC, Boca Raton, FL. ISBN 978-0-41204231-2.

Simon J. Sheather and M. Chris Jones (1991) A reliable data-based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society, Series B (Methodological), 53(3), 683--690.

Matt P. Wand and M. Chris Jones (1995) Kernel Smoothing. Chapman and Hall, London, UK. ISBN 978-0-412-55270-1.}

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
bootkern2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
           expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
