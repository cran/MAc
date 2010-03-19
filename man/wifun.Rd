\name{Wifun}
\title{Weights added to Meta Data}
\alias{Wifun}
\description{Adds weights to the meta-analysis data set.}
\usage{
Wifun(meta)
}
\arguments{
\item{meta}{\code{data.frame} with id, n (sample size), and r (correlation coefficient) for each study. The function requires that the column name for these values be exactly as stated here.
}
}
\value{Adds fixed and random-effects weights and confidence intervals to meta data set.
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{data}