\name{CatCompr}
\alias{CatCompr}
\title{Direct Categorical Moderator Comparison (Random Effects)
}
\description{Function for a planned comparison between two levels of a moderator under a random effects model.
}
\usage{
CatCompr(meta, mod, x1, x2, method = "post.hoc")
}
\arguments{
  \item{meta}{data.frame with correlation coefficients and sample size for each study.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
  \item{x1}{One level of categorical moderator.
}
  \item{x2}{Comparison level of same categorical moderator.
}
  \item{method}{post.hoc assumes the comparison was not planned prior to conducting the meta analysis. A planned comparison assumes the researcher has planned to compare these levels of the categorical moderator a priori. Default is post.hoc. 
}
}
\details{See page 280-288 in Cooper et al. (2009) for the computations used in this function.
}
\value{
\item{Mean1}{Mean effect size of the first level of the moderator. compared.
}
\item{Mean2}{Mean effect size of the second level of the moderator. compared.
}
\item{d}{ Mean difference between the two groups divided by its standard deviation.
}
\item{Var d}{Variance of d.
}
\item{LowerLimits}{ Lower 95\% confidence interval.
}
\item{UpperLimits}{ Upper 95\% confidence interval.
}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta analysis (2nd edition). New York Russell Sage Foundation.
} 
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatCompf}},
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

CatCompr(data,data$mod1,1,3) # compares level 1 vs level 3 of moderator. 
}
\keyword{models}
\keyword{regression}


