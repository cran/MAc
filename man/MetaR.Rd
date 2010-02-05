\name{MetaR}
\alias{MetaR}
\title{Meta-Analysis Aggregation
}
\description{
Computes within-study aggregation (see Hunter & Schmidt, 2004) and adds fixed and random effects weights for each study in the meta-analysis.
}
\usage{
MetaR(meta)
}
\arguments{
  \item{meta}{data.frame with r (correlation coefficients) and n (sample size) for each study.}
}
\value{Outputs a data.frame with aggregated effect sizes where each study is reduced to one row per study by a weighted average formula. This formula is based on Hunter and Schmidts (2004) approach to aggregation of dependent effect sizes for correlations (see chapter 10, pp. 435-8). Further, this function will add all necessary weights to each study in order to conduct all appropriate analyses. These weights are derived from a statistical transformation of each r (correlation coefficient) to a Fishers z (see chapter 14, Cooper et al., 2009).
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
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
data<-data.frame(id,n,r,mod1)

# Example

MetaR(data)   
}
\keyword{ models }

