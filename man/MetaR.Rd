\name{MetaR}
\alias{MetaR}
\title{Meta-Analysis Aggregation
}
\description{
Computes within-study aggregation (see Hunter & Schmidt, 2004) and adds fixed and random effects weights for each study in the meta-analysis.
}
\usage{
MetaR(meta, cor = .50)
}
\arguments{
  \item{meta}{\code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
}
\value{Outputs a \code{data.frame} with aggregated effect sizes where each study is reduced to one row per study by a weighted average formula. This formula is based on Hunter and Schmidt's (2004) approach to aggregation of dependent effect sizes for correlations (see chapter 10, pp. 435-8). Further, this function will add all necessary weights to each study in order to conduct all appropriate analyses. These weights are derived from a statistical transformation of each r (correlation coefficient) to a Fisher's z (Shadish & Haddock, 2009).
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 

Hunter, J. E., Schmidt, F. L. (2004). \emph{Methods of meta-analysis} (2nd edition). Thousand Oaks, CA: Sage. 

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.

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

