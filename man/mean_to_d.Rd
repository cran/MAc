\name{mean_to_d}
\alias{mean_to_d}
\title{Means to Standardized Mean Difference
}
\description{Converts raw mean scores reported in the primary study to a standardized mean difference (d).
}
\usage{
mean_to_d(m.1, m.2, sd.1, sd.2, n.1, n.2)
}
\arguments{
  \item{m.1}{Mean of group one.
}
  \item{m.2}{Mean of group two.
}
  \item{sd.1}{Standard deviation of group one.
}
  \item{sd.2}{Standard deviation of group two.
}
  \item{n.1}{Sample size of group one. 
}
  \item{n.2}{Sample size of group two. 
}
}
\value{
\item{d}{Standardized mean difference (d).}
\item{var.d }{Variance of d.}

}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
}
\seealso{
\code{\link{d_to_g}},
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
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\keyword{ arith }