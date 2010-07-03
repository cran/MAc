\name{wgts}
\title{Weights added to Meta Data}
\alias{wgts}
\description{Adds weights to the meta-analysis data set.}
\usage{
wgts(es, var.es, data)
}
\arguments{
\item{es}{r or z' effect size. 
}
\item{var.es}{Variance of es 
}
\item{data}{\code{data.frame} with values above. 
}
}
\value{Adds fixed and random-effects weights and confidence intervals to meta data set.
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{data}