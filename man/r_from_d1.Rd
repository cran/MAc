\name{r_from_d1}
\alias{r_from_d1}
\title{Correlation from Mean Difference II
}
\description{Converts d (mean difference) reported in primary study to r (correlation coefficient) where sample size of the treatment group  is not equal to the sample size of the comparison group (see section 12.5.4, Cooper et al., 2009; pp. 234)
}
\usage{
r_from_d1(d, n.1, n.2, var.d)
}
\arguments{
  \item{d}{Mean difference statistic (d) reported in primary study.
}
  \item{n.1}{sample size of group 1 reported in primary study.
}
  \item{n.2}{sample size of group 2 reported in primary study.
}
  \item{var.d}{variance of d reported in primary study.
}
}
\value{ Computes a correlation coefficient (r)}
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

