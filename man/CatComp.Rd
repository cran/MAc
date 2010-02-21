\name{CatComp}
\alias{CatComp}
\title{Direct Categorical Moderator Comparison 
}
\description{Function for a planned comparison between two levels of a moderator under a fixed and random effects model.
}
\usage{
CatComp(meta, mod, x1, x2, method = "post.hoc1")
}
\arguments{
  \item{meta}{\code{data.frame} with correlation coefficients and sample size for each study. The function requires that the column name for the correlation coefficient is 'r' and the sample size is 'n'.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
}
  \item{x1}{One level of categorical moderator.
}
  \item{x2}{Comparison level of same categorical moderator.
}
  \item{method}{\code{post.hoc1} assumes the comparison was not planned prior to conducting the meta analysis. The other option, \code{planned}, assumes the researcher planned to conduct the analysis a priori. Default is \code{post.hoc1} using Tukey's Honest Significant Differences method. \code{post.hoc2} uses the Scheffe post hoc statistical method. 
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{
\item{d}{ Mean difference between the two groups divided by its standard deviation (Standardized mean difference).
}
\item{Var(d)}{Variance of d (standardized mean difference).
}
\item{CI.lower}{ Lower 95\% confidence interval.
}
\item{CI.lower}{ Upper 95\% confidence interval.
}
}
\references{ Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation. 
} 
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{CatMod}},
\code{\link{CatModGraph}}
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

CatComp(datas,datas$mod1,1,3) # compares level 1 vs level 3 of moderator. 
}
\keyword{models}
\keyword{regression}
