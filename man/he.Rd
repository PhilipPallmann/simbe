\name{he}

\alias{he}

\title{Simultaneous confidence bounds for a multivariate normal mean using He's empirical Bayes method}

\description{Computation of simultaneous confidence limits around multivariate normal mean estimates (with shrinkage towards 0 for more than 2 parameters) using the empirical Bayes method of He (1992).}

\usage{he(dat, alpha=0.1)}

\arguments{
  \item{dat}{A data frame with as many rows as independent observations and as many columns as endpoints.}
  \item{alpha}{A numerical value specifying the type I error level.}
}

\details{truncations etc. The version implemented here is given as equation (3) in Qiu & Hwang (2007).}

\value{A matrix of the parameter estimates and simultaneous confidence intervals.}

\references{Kun He (1992) Parametric empirical Bayes confidence intervals based on James-Stein estimator. Statistics & Decisions, 10(1--2), 121--132.

Jing Qiu and J. T. Gene Hwang (2007) Sharp simultaneous confidence intervals for the means of selected populations with application to microarray data analysis. Biometrics, 63(3), 767--776.}

\author{Philip Pallmann}

\note{xxx}

\section{Warning}{Functionality for 3 or more parameters has not yet been thoroughly tested.}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 3), sigma=diag(3) * 0.02)
colnames(mydata) <- c("AUC", "Cmax", "xyz")
he(mydata)

# compare with standard CIs
}

\keyword{multivariate}
