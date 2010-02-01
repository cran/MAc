\name{MAregGraph}
\alias{MAregGraph}
\title{Meta Regression Scatterplot
}
\description{Outputs a scatterplot from a fixed or random effects meta regression (continuous and/or categorical). 
}
\usage{
MAregGraph(meta, mod, method = "random", modname = NULL, title = NULL, ylim = c(0, 1))
}
\arguments{
  \item{meta}{ data.frame with r (correlation coefficients) and n (sample size) for each study.
}
  \item{mod}{ Moderator variable used for meta regression.  
}
  \item{method}{ Model used, either random or fixed effects. Default is random.
}
  \item{modname}{Name of moderator to appear on x axis of plot. Default is NULL.
}
  \item{title}{Plot title. Default is NULL.
}
  \item{ylim}{Limits of y axis with the first arguments minimal value and second maximum value. Default is c(0, 1).
}
}
\value{ Scatterplot with fixed or random effects regression line where size of points are based on study weights, where the more precise studies are larger. The ggplot2 package outputs the rich graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{MAreg1}},
\code{\link{MAreg2}}
}
\examples{
# Data Set III (Cooper et al., 2009; pp. 265 & 587)

id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
data<-data.frame(id,n,r,mod1)

# Example

MAregGraph(data, data$mod1)   
}
\keyword{ aplot }

