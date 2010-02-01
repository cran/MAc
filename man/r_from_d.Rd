\name{r_from_d}
\alias{r_from_d}
\title{Correlation from Mean Difference
}
\description{Converts d (mean difference) to r where n.tmt = n.comparison (see section 12.5.4, Cooper, et al., 2009; pp. 234)
}
\usage{
r_from_d(d, var.d, a = 4)
}
\arguments{
  \item{d}{Mean difference statistic (d) reported in primary study
 }
  \item{var.d}{Variance of d reported in primary study
}
  \item{a}{Used to compute correlation. Default is 4 (and should not be changed).
}
}
\value{Outputs a correlation coefficient (r)}
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

