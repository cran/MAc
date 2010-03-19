\name{agg_r}
\alias{agg_r}
\title{Intermediate level function}
\description{Intermediate level function to aggregate data and create new data.frame. Required inputs are data.frame with r (correlation coefficients) and n (sample size).
}
\usage{
agg_r(meta, cor = .50)
}
\arguments{
  \item{meta}{data.frame with r (correlation coefficients) and n (sample size)for each study.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
}
\value{Aggregated correlations
}
\references{Hunter, J. E., Schmidt, F. L. (2004). Methods of meta-analysis (2nd edition). Thousand Oaks, CA: Sage.
}
\author{AC Del Re & William T. Hoyt
}
\keyword{internal}



