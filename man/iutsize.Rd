\name{iutsize}
\alias{iutsize}

\title{Actual size of the x}

\description{intersection union test}

\usage{iutsize(p, n, alpha=0.1, sim=1e5)}

\arguments{
  \item{p}{xxx}
  \item{n}{xxx}
  \item{alpha}{xxx}
  \item{sim}{xxx}
}

\details{}

\value{A numeric value giving the actual size of the test.}

\references{Roger L. Berger (1982) Multiparameter hypothesis testing and acceptance sampling. Technometrics, 24(4), 295--300.

Weizhen Wang, J. T. Gene Hwang & Anirban DasGupta (1999) Statistical tests for multivariate bioequivalence. Biometrika, 86(2), 395--402.}

\author{Philip Pallmann}

\note{xxx}

\seealso{x}

\examples{
# For p=1 we get the well-known result that the
# ordinary 90% CI corresponds to the TOST at 5%:

iutsize(p=1, n=20)

# With increasing dimension the test gets conservative:

iutsize(p=2, n=20)
iutsize(p=3, n=20)

# The conservativeness also depends on sample size:

iutsize(p=2, n=10)
iutsize(p=2, n=1000)
}

\keyword{ ~kwd1 }
