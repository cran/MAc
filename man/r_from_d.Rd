\name{r_from_d}
\alias{r_from_d}
\title{Correlation from Mean Difference
}
\description{Converts d (mean difference) to r where n.tmt = n.comparison (see section 12.5.4, Borenstein, 2009; pp. 234).
}
\usage{
r_from_d(d, var.d, a = 4)
}
\arguments{
  \item{d}{Mean difference statistic (d) reported in primary study
 }
  \item{var.d}{Variance of d reported in primary study
}
  \item{a}{Used to compute correlation. Default is 4.
}
}
\value{Outputs a correlation coefficient (r)
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{ 
\code{\link{r_from_d}},
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\keyword{ arith }

