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

The outcome variable used in the initial regression should be a Fishers z (which will be coverted back to r in the final output). The predictors (moderator(s)) can be either continuous or categorical. Weight the regression by either the fixed or random effect weight (e.g., fixed equals \code{data$wi} and random equals \code{data$wi.tau}). See \code{ComplData} function to assist in reducing the dataset, if applicable. See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
}
\value{
\item{b}{ Unstandardized regression coefficient.
} 
\item{SE}{ Standard error of b coefficient (adjusted SE value).
}
\item{t}{ t-value (adjusted value).
}
\item{CI.Lower}{ Lower 95\% confidence interval.
}
\item{CI.Upper}{ Upper 95\% confidence interval.
}
\item{p}{ Significance level.
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

