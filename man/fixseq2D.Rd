\name{fixseq2D}

\alias{fixseq2D}

\title{Plot bivariate (simultaneous) confidence intervals corresponding to two one-sided tests for two parameters}

\description{xxx}

\usage{fixseq2D(dat, alpha=0.1, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Fixed Sequence")}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
  \item{equi}{xxx}
  \item{plotrange}{xxx}
  \item{axnames}{xxx}
  \item{main}{xxx}
}

\details{xxx}

\value{xxx}

\references{Willi Maurer, Ludwig A. Hothorn, and Walter Lehmacher (1995) Multiple comparisons in drug clinical trials and preclinical assays: a priori ordered hypotheses. In: Joachim Vollmar (editor), Biometrie in der Chemisch-Pharmazeutischen Industrie, vol. 6, pp. 3--18. Fischer-Verlag, Stuttgart, Germany.

Strassburger Bretz

Guilbaud}

\author{Philip Pallmann}

\note{xxx}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=diag(2) * 0.05)
fixseq2D(mydata, axnames=c(expression(paste(Delta, " log(AUC)")),
         expression(paste(Delta, " log(Cmax)"))))

# also situations where the other cases happen
}

\keyword{multivariate}
