\name{hotelling2D}

\alias{hotelling2D}

\title{Plot bivariate simultaneous confidence region of Wang et al. (1999)}

\description{xxx}

\usage{hotelling2D(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3),
            axnames=NULL, main="Hotelling")}

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

\references{Harold Hotelling (1931) The generalization of Student's ratio. The Annals of Mathematical Statistics, 2(3), 360--378.

Weizhen Wang, J. T. Gene Hwang & Anirban DasGupta (1999) Statistical tests for multivariate bioequivalence. Biometrika, 86(2), 395--402.}

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
hotelling2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
          expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
