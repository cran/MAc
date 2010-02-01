\name{cat_sum1R}
\alias{cat_sum1R}
\title{Intermediate level function}
\description{Intermediate level function for categorical moderation. sums relevant es information by each moderator level. 
}
\usage{
cat_sum1R(meta, mod)
}
\arguments{
  \item{meta}{ data.frame with r (correlation coefficients) and n (sample size)for each study.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
}
\value{Random effects means and weights for each level of moderator. 
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta-analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{ AC Del Re & William T. Hoyt
}
\keyword{internal}

