\name{ComplData}
\alias{ComplData}
\title{Complete Dataset
}
\description{ Function to reduce data set with complete data for 1 predictor. Will accept a \code{data.frame} with multiple rows per study and reduce the \code{data.frame} to one effect size per study based on recommended aggregation procedures for correlations (Hunter & Schmidt, 2004; pp. 435-8).
}
\usage{
ComplData(meta, mod, type = "independent", cor = .50)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study. The function requires that the column name for the correlation coefficient is 'r' and the sample size is 'n'.
}
  \item{mod}{ Moderator variable wanting to be kept for further analysis. 
}
  \item{type}{ Type of data aggregation. Options are \code{independent} or \code{dependent}. Choose \code{independent} if a fully independent data set is desired, for which the function will randomly select one factor level of the within-study moderator (if there are > 1 coded from the primary studies, e.g., the primary study reports early, mid, and late symptom change from CBT treatment). Choose \code{dependent} if all factor levels of the within-study moderator are to be kept in the data set. Default is \code{independent}.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. recommended procedures (1997).
}
}
\value{ Reduced \code{data.frame} (with complete data) for the moderator entered into the function. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Raudenbush (2009). Analyzing effect sizes: Random-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 295-316). New York: Russell Sage Foundation.  

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.
}
\examples{
# Sample data from Data Set III (Raudenbush, 2009; pp. 300)

id<-1:20
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
datas<-data.frame(id,n,r,mod1) 

# Example

ComplData(datas, datas$mod1)
}
\keyword{data}

