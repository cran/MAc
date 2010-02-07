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
  \item{meta}{\code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{var}{\code{weighted} or \code{unweighted}. Default is \code{weighted}. Use the \code{unweighted} variance method only if Q is rejected and is very large relative to the number of studies in the meta-analysis. 
}
}
\value{
Fixed and random effects:
 
\item{K}{ Number of studies in the meta-analysis.
}
\item{ES}{ Overall Omnibus effect size for the meta-analysis.
}
\item{Var(ES)}{ Variance of Omnibus effect size.
}
\item{StdError}{ Square root of variance.
}
\item{LowerLimit}{ Lower 95\% confidence interval.
}
\item{UpperLimit}{ Upper 95\% confidence interval.
}
\item{Z-value}{ Z-score (standardized value).
}
\item{P-value}{ Significance level.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{Q(df)}{ Degrees of freedom for Q-statistic.
}
\item{P-hetero}{ Q-statistic p-value (assesses overall heterogeneity between studies).
}
\item{I-squared}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
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

