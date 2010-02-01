\name{r_from_t}
\alias{r_from_t}
\title{Correlation from T-test
}
\description{Compute correlation from t-value reported in the primary study}
\usage{
r_from_t(t, n)
}
\arguments{
  \item{t}{ t-statistic value reported in primary study.
}
  \item{n}{ Sample size reported in primary study.
}
}
\value{Outputs the correlation coefficient (r).
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta-analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{ 
\code{\link{r_from_d}},
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\keyword{ arith }

