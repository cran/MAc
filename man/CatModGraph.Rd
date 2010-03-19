\name{CatModGraph}
\alias{CatModGraph}
\title{Categorical Moderator Graph
}
\description{ Outputs a rich and detailed boxplot graphic for each level of the specified moderator (under a fixed or random effects model). 
}
\usage{
CatModGraph(meta, mod, method = "random", modname = NULL,title = NULL)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study. The function requires that the column name for the correlation coefficient is 'r' and the sample size is 'n'.
}
  \item{mod}{ Categorical moderator variable used for analysis. 
}
  \item{method}{ Model used, either \code{random} or \code{fixed} effects model. Default is \code{random}. 
}
  \item{modname}{ Name of moderator to appear on x-axis of plot. Default is NULL. 
}
  \item{title}{ Plot title. Default is NULL. 
}
}
\value{ Boxplot graph with median, inter-quartile range, max, min, and outliers from a fixed or random effects categorical moderator analysis. Places jitter points (for each study) on the boxplots. The size of each point (representing a study in the analysis) are based on study weights where more precise studies are larger. The ggplot2 package outputs the rich graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatComp}},
\code{\link{CatMod}},
}
\examples{
# Data Set III (Cooper et al., 2009; pp. 265 & 587)

id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<- factor(rep(1:4,5))
datas<-data.frame(id,n,r,mod1)

# Example

\dontrun{CatModGraph(datas,datas$mod1)}
}   
\keyword{aplot}

