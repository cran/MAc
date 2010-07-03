\name{d_to_g}
\alias{d_to_g}
\title{
Standardized Mean Difference (d) Statistic to Unbiased Standardized Mean Difference (g)
}
\description{
Eliminates the small upward bias of d to provide an unbiased estimate of the population effect size parameter (g). This procedure will compute g for a single value of d. 
}
\usage{
d_to_g(d, var.d, n.1, n.2)
}
\arguments{
  \item{d}{Standardized mean difference statistic (d has a slight bias).
}
  \item{var.d}{Variance of d.
}
  \item{n.1}{Sample size of treatment group.
}
  \item{n.2}{Sample size of comparison group.
}
}
\value{
 \item{g }{Unbiased Standardized mean difference statistic. }
 \item{var.g }{Variance of g.}
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}

\keyword{ arith }