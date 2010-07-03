\name{r_from_chi}
\alias{r_from_chi}
\title{Chi-Squared to Correlation
}
\description{Converting Chi-squared statistic (with 1 degree of freedom) reported in primary study to r
}
\usage{
r_from_chi(chi.sq, n)
}
\arguments{
  \item{chi.sq}{Chi squared statistic from primary study.
}
  \item{n}{Sample size in primary study.
}
}
\value{Computes correlation coefficient (r).
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
\code{\link{r_from_d}},
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\keyword{ arith }

