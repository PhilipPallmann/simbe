\name{tseng2D}

\alias{tseng2D}

\title{Plot bivariate simultaneous confidence region of Tseng (2002)}

\description{xxx}

\usage{tseng2D(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
}

\details{xxx}

\value{xxx
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{Yu-Ling Tseng (2002) Optimal confidence sets for testing average bioequivalence. Test, 11(1), 127--141.}

\author{Philip Pallmann}

\note{When the estimates are close to zero, the confidence region may not exist.}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=diag(2) * 0.05)
tseng2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
        expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
