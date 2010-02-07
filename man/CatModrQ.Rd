\name{CatModrQ}
\alias{CatModrQ}
\title{Categorical Moderator Q-statistic (Random Effect) 
}
\description{Q test for homogeneity between levels of moderator under a random effects model.
}
\usage{
CatModrQ(meta, mod)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{mod}{ Categorical moderator variable used for moderator analysis.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{ 
Random effects moderator:

\item{Q}{ Q-statistic (measure of homogeneity).
} 
\item{Qw}{ Q-within (measure of within-group heterogeneity).
}
\item{df(Qw)}{ Degrees of freedom for Q-within.
}
\item{P-value(within)}{ Q-within p-value (for homogeneity)
}
\item{Qb}{ Q-between (measure of between-group heterogeneity)
}
\item{df(Qb)}{ Degrees of freedom for Q-between.
}
\item{P-value(between)}{ Q-between p-value (for homogeneity)
} 
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatModr}},
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

CatModrQ(data,data$mod1)
}
\keyword{models}

