\name{CatModr}
\alias{CatModr}
\title{Categorical Moderator Analysis (Random Effects) 
}
\description{Computes single predictor categorical moderator analysis.
}
\usage{
CatModr(meta, mod)
}
\arguments{
  \item{meta}{data.frame with r (correlation coefficients) and n (sample size) for each study.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
}
\value{
Random effects moderator:

\item{Mod}{ Level of the categorical moderator.
} 
\item{K}{ Number of studies for each level of the moderator.
}
\item{ES}{ Mean effect size of each level of the moderator.
}
\item{LowerLimit}{ Lower 95\% confidence interval.
}
\item{UpperLimit}{ Upper 95\% confidence interval.
}
\item{Z-value}{ Z-score (standardized value).
}
\item{P-value}{ Significance level.
}
\item{Var}{ Variance of effect size.
}
\item{StdError}{ Sqrt of variance.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{Q(df)}{ Degrees of freedom for Q-statistic.
}
\item{I-squared}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Cooper, et al., 2009; pp. 263).
}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta-analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatCompr}},
\code{\link{CatModf}},
\code{\link{CatModrQ}},
\code{\link{CatModfQ}},
\code{\link{CatModGraph}}
}
\examples{
# Data Set III (Cooper et al., 2009; pp. 265 & 587)

id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<-factor(c(rep(c(1,2,3,4),5)))   
data<-data.frame(id,n,r,mod1)

# Example

CatModr(data,data$mod1)   
}
\keyword{models}

