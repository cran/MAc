\name{MRfit}
\alias{MRfit}
\title{Meta-Regression Model Fit 
}
\description{ Computes the model fit of a moderator analysis comparing an empty model to one with x number of predictors.
}
\usage{
MRfit(...)
}
\arguments{
  \item{...}{ Any number of weighted linear regressions (saved as objects, see details in example section below).
}
}
\details{ Running this function requires the researcher to conduct a weighted hierarchical linear regression while saving each block as an object (see example below).
The outcome variable is a Fishers z and the predictor moderator(s) can be either continuous or categorical. Weight the regression by either the fixed or random effect weight (e.g., fixed equals \code{data$wi} and random equals \code{data$wi.tau}). See \code{ComplData} function to assist in reducing the dataset, if applicable. See Konstantopoulos & Hedges (2009; pp. 289-290) for the computations used in this function.
}
\value{
\item{df.Q}{ Q degrees of freedom.
} 
\item{Qe}{ Q-error, weighted residual sum of squares.
}
\item{predictors}{ Number of predictors in the model for each block.
}
\item{Q}{ Q-statistic, indicating overall heterogeneity.
}
\item{p}{ Significance level.
}
\item{R^2}{ R-squared, the coefficient of determination, is the proportion of variability in a data set that is accounted for by the statistical model.
}
\item{R^2.change}{ Change in R-squared from one model to the next (from one block to the next).

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

reg0 <- lm(data$z ~ 1, weights= data$wi.tau)  # empty model
reg1 <- lm(data$z ~ data$mod2, weights= data$wi.tau)  # 1 predictor model
reg2 <- lm(data$z ~ data$mod2 + data$mod1, weights= data$wi.tau)  # 2 predictor model
MRfit(reg0, reg1, reg2) 
}  
\keyword{ model }
\keyword{ regression }

