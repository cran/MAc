\name{r_from_d1}
\alias{r_from_d1}
\title{Correlation from Mean Difference II
}
\description{Converts d (mean difference) reported in primary study to r (correlation coefficient) where sample size of the treatment group  is not equal to the sample size of the comparison group (see section 12.5.4, Borenstein, 2009; pp. 234).
}
\usage{
r_from_d1(d, n.1, n.2, var.d)
}
\arguments{
  \item{d}{Mean difference statistic (d) reported in primary study.
}
  \item{n.1}{sample size of group one reported in primary study.
}
  \item{n.2}{sample size of group two reported in primary study.
}
  \item{var.d}{variance of d reported in primary study.
}
}
\value{ Computes a correlation coefficient (r).}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{d_to_g}},
\code{\link{mean_to_d}},
\code{\link{mean_to_d2}},
\code{\link{t_to_d}},
\code{\link{f_to_d}},
\code{\link{p_to_d1}},
\code{\link{p_to_d2}}, 
\code{\link{ancova_to_d1}}, 
\code{\link{ancova_to_d2}}, 
\code{\link{tt.ancova_to_d}}, 
\code{\link{f.ancova_to_d}}, 
\code{\link{r_to_d}}, 
\code{\link{p.ancova_to_d1}},
\code{\link{p.ancova_to_d2}},
\code{\link{lor_to_d}},
\code{\link{prop_to_or}}, 
\code{\link{prop_to_d}}, 
\code{\link{r_from_chi}}, 
\code{\link{r_from_d}},
\code{\link{r_from_t}}
}
\keyword{ arith }

