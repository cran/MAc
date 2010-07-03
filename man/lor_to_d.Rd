\name{lor_to_d}
\alias{lor_to_d}
\title{Log Odds Ratio to Standardized Mean Difference (d)
}
\description{Converts a log odds ratio reported in the primary study to a standardized mean difference (d).
}
\usage{
lor_to_d(lor, var.lor)
}
\arguments{
  \item{lor}{Log odds ratio reported in the primary study.
}
  \item{var.lor}{Variance of the log odds ratio.
}
}
\value{
\item{d}{Standardized mean difference (d).
}
\item{var_d }{Variance of d.
}

}
\author{ AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Borenstein (2009). Effect sizes for continuous data. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.
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
\code{\link{prop_to_or}}, 
\code{\link{prop_to_d}}, 
\code{\link{r_from_chi}}, 
\code{\link{r_from_d}},
\code{\link{r_from_d1}},
\code{\link{r_from_t}}
}
\keyword{ arith }