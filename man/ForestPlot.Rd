\name{ForestPlot}
\alias{ForestPlot}
\title{Forest Plot
}
\description{Forest plot for meta-analysis}
\usage{
ForestPlot(meta, method = "random", title = NULL)
}
\arguments{
  \item{meta}{data.frame with r (correlation coefficients) and n (sample size) for each study.
}
  \item{method}{Model used, either "random" or "fixed" effects. Default is "random".
}
  \item{title}{Plot title. Default is NULL.
}
}
\value{Forest plot with omnibus effect size (fixed or random), point for each study where size of point is based on the study's precision (based primarily on sample size) and 95\% confidence intervals. The ggplot2 package outputs the rich graphics.  
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta-analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
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

ForestPlot(data)   
}
\keyword{aplot}

