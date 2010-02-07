\name{r_from_t}
\alias{r_from_t}
\title{Correlation from t-test
}
\description{Compute correlation from t-value reported in the primary study (see section 12.5.4, Borenstein, 2009; pp. 234).
}
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
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
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

