\name{fixseq}

\alias{fixseq}

\title{Simultaneous confidence bounds for a fixed sequence procedure}

\description{only 2D
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
fixseq(dat, alpha=0.1)
}

\arguments{
  \item{dat}{xxx}
  \item{alpha}{xxx}
}

\details{
xxx
}

\value{A matrix of the parameter estimates and (marginalized) simultaneous confidence intervals.}

\references{
Willi Maurer, Ludwig A. Hothorn, and Walter Lehmacher (1995) Multiple comparisons in drug clinical trials and preclinical assays: a priori ordered hypotheses. In: Joachim Vollmar (editor), Biometrie in der Chemisch-Pharmazeutischen Industrie, vol. 6, pp. 3--18. Fischer-Verlag, Stuttgart, Germany.

Strassburger Bretz

Guilbaud
}

\author{Philip Pallmann}

\note{xxx not for p greater than 2
%%  ~~further notes~~
}

\seealso{xxx
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
library("mvtnorm")
mydata <- rmvnorm(50, mean=rep(0.1, 2), sigma=diag(2) * 0.05)
colnames(mydata) <- c("AUC", "Cmax")
fixseq(mydata)

# also situations where the other cases happen
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{ ~kwd1 }
