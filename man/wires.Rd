\name{wires}

\alias{wires}

\docType{data}

\title{
Quality control study of wire connections
}

\description{
Data from a study in quality control assessing the breaking strengths of 20 wire connections between a semiconductor wafer and a terminal post (King 1971).
}

\usage{data("wires")}

\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{Strength}}{A numeric vector giving the strength at which the connection failed.}
    \item{\code{Failure}}{A factor with levels \code{b} and \code{w} specifying whether the bond or wire failed.}
  }
}

\details{
The data were taken from Table 4.1 of Nelson (2004).
}

\source{
King 1971
}

\references{
Wayne B. Nelson (2004) Applied Life Data Analysis. Wiley, Hoboken, NJ.
}

\examples{
data(wires)

# Simultaneous confidence regions for the mean and variance
plotMV2D(wires$Strength, method="mood", main="Mood")
plotMV2D(wires$Strength, method="lrt", main="LRT")
}

\keyword{datasets}
