\name{FunnelPlot}
\alias{FunnelPlot}
\title{Funnel Plot
}
\description{Outputs a funnel plot from a fixed or random effects omnibus analysis to visually assess meta-analytic data for publication bias.
}
\usage{
FunnelPlot(meta, method = "random", title = NULL)
}
\arguments{
  \item{meta}{\code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{method}{Plot Model used is either \code{random} or \code{fixed} effects. Default is \code{random}.
}
  \item{title}{Plot title. Default is NULL.
}
}
\value{Outputs a funnel plot with omnibus effect size (fixed or random), visual point for each study, where size of point is based on the study's precision (based primarily on sample size), and predicted standard error lines drawn to assist in assessing for publication bias in the meta-analysis. The ggplot2 package outputs the rich graphics. 
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 
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
datas<-data.frame(id, n, r, mod1)

# Example

\dontrun{FunnelPlot(datas)}
}
\keyword{aplot}
