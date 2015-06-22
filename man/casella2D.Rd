\name{casella2D}

\alias{casella2D}

\title{Plot bivariate simultaneous confidence region of Casella & Hwang (1983)}

\description{xxx}

\usage{casella2D(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL)}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{steps}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
}

\details{xxx}

\value{xxx}

\references{George Casella & Jiunn Tzon Hwang (1983) Empirical Bayes confidence sets for the mean of a multivariate normal distribution. Journal of the American Statistical Association, 78(383), 688--698.}

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
casella2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
          expression(paste(Delta, " log(Cmax)"))))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
