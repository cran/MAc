\name{MAreg1}
\alias{MAreg1}
\title{ Meta-Regression (Single Predictor Model)
}
\description{ Meta-regression function for a single predictor model. This function corrects for inaccuracies in the standard errors produced from an OLS regression model (Konstantopoulos & Hedges, 2009; pp. 289-290).
}
\usage{
MAreg1(meta, mod, method = "random")
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{mod}{ Moderator variable used for meta-regression.  
}
  \item{method}{ Model used, either \code{random} or \code{fixed effects}. Default is \code{random}.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
}
\value{
\item{Beta}{ Beta coefficient.
} 
\item{StdError}{ Standard Error of beta coefficient (adjusted StdError value).
}
\item{Tvalue}{ T-value (adjusted value).
}
\item{LowerLimit}{ Lower 95\% confidence interval.
}
\item{UpperLimit}{ Upper 95\% confidence interval.
}
\item{P-value}{ Significance level.
}
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{MAreg1}},
\code{\link{MAreg2}},
\code{\link{MAregGraph}}
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

MAreg1(data, data$mod1)   
}
\keyword{ models }
\keyword{ regression }


