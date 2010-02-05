\name{MAreg2}
\alias{MAreg2}
\title{Meta-Regression (Multiple Predictor Model) 
}
\description{ Computes multiple predictor fixed or random effects meta-regression (continuous and/or categorical).
}
\usage{
MAreg2(reg)
}
\arguments{
  \item{reg}{ Weighted linear regression (saved as an object, see details below).
}
}
\details{ Running this function requires the researcher to conduct a weighted linear regression while saving it as an object, e.g.:
  
\code{reg <- lm(data$z ~ data$mod8 + data$mod1, weights= data$wi.tau)}

The outcome variable is a Fishers z and the predictor moderator(s) can be either continuous or categorical. Weight the regression by either the fixed or random effect weight (e.g., fixed equals \code{data$wi} and random equals \code{data$wi.tau}). See \code{ComplData} function to assist in reducing the dataset, if applicable. See page 289 and 290 in Cooper et al. (2009) for the computations used in this function.
}
\value{
\item{Beta}{ Beta coefficient.
} 
\item{StdError}{ Standard Error of coefficient (adjusted value).
}
\item{Tvalue}{ T value (adjusted value).
}
\item{LowerLimit}{ Lower 95\% confidence interval.
}
\item{UpperLimit}{ Upper 95\% confidence interval.
}
\item{P-value}{ Significance level.
}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta analysis (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{MAreg1}},
\code{\link{MAreg2}},
\code{\link{MAregGraph}},
\code{\link{ComplData}}

}
\examples{
# Data Set III (Cooper et al., 2009; pp. 265 & 587)

require(MAc)
id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
mod2<-factor(c(rep(c(1,2,3,4),5))) 
data<-data.frame(id, n, r, mod1, mod2)
data$z <- r_to_z(data$r) # r to z function in this package

# Example

regr <- lm(data$z ~ data$mod1 + data$mod2, weights= data$wi)
 
MAreg2(regr) # Function will correct standard errors
}  
\keyword{ model }
\keyword{ regression }

