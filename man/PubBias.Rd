\name{PubBias}
\alias{PubBias}
\title{Assess Publication Bias
}
\description{Assess for publication bias in the meta-analytic data}
\usage{
PubBias(meta)
}
\arguments{
 \item{meta}{\code{data.frame} with id and z.score (standardized z-value). The function requires that the column name for these values be exactly as stated here.
}
}
\value{
 \item{k}{Number of studies.}
 \item{Z}{Overall z-value for data set.}
 \item{k0}{Number of studies needed to include with effect size = 0 (null) in order for the p > .05 (null hypothesis retained).}
  \item{k.per}{Number of missing studies for every observed study for the overall effect to be nullified.}
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Sutton, A. J. (2009). Publication bias. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 435-452). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{FunnelPlot}},
}
\keyword{ data }