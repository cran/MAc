\name{OmnibusES}
\alias{OmnibusES}
\title{Omnibus Effect Size (Fixed and Random Effects) 
}
\description{Computes fixed and random effects omnibus effect size for correlations. 
}
\usage{
OmnibusES(meta, var = "weighted")
}
\arguments{
  \item{meta}{\code{data.frame} with r (correlation coefficients) and n (sample size) for each study. The function requires that the column name for the correlation coefficient is 'r' and the sample size is 'n'.
}
  \item{var}{\code{weighted} or \code{unweighted}. Default is \code{weighted}. Use the \code{unweighted} variance method only if Q is rejected and is very large relative to the number of studies in the meta-analysis. 
}
}
\value{
Fixed and random effects:
 
\item{k}{ Number of studies in the meta-analysis.
}
\item{ES}{ Overall Omnibus effect size for the meta-analysis.
}
\item{var.ES}{ Variance of Omnibus effect size.
}
\item{SE}{ Square root of variance.
}
\item{CI.lower}{ Lower 95\% confidence interval.
}
\item{CI.upper}{ Upper 95\% confidence interval.
}
\item{Z}{ Z-score (standardized value).
}
\item{p}{ Significance level.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df}{ Degrees of freedom for Q-statistic.
}
\item{p.h}{ Q-statistic p-value (assesses overall homogeneity between studies).
}
\item{I2}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
}
}
\references{ Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
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

OmnibusES(data)
}
\keyword{ models }

