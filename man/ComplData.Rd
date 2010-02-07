\name{ComplData}
\alias{ComplData}
\title{Complete Dataset (by Moderator(s))
}
\description{ Function to reduce data set with complete data for \emph{x} number of predictors in a multivariate moderator analysis. Will accept a \code{data.frame} with multiple rows per study and reduce the \code{data.frame} to one effect size per study based on recommended aggregation procedures for correlations (Hunter & Schmidt, 2004; pp. 435-8).
}
\usage{
ComplData(meta, mod1, mod2 = NULL, mod3 = NULL, mod4 = NULL,
 mod5 = NULL, predictors = 1)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{mod1}{ Moderator variable wanting to be kept for further analysis. 
}
  \item{mod2}{ Moderator variable wanting to be kept for further analysis. Default is NULL.
}
  \item{mod3}{ Moderator variable wanting to be kept for further analysis. Default is NULL.
}
  \item{mod4}{ Moderator variable wanting to be kept for further analysis. Default is NULL.
}
  \item{mod5}{ Moderator variable wanting to be kept for further analysis. Default is NULL.
}
  \item{predictors}{ Number of moderators specified in function (to keep in reduced \code{data.frame}).
}
}
\value{ Reduced \code{data.frame} (with complete data) for each moderator entered into the function. This is primarily used as a convenience function for conducting multivariate meta-regressions or when desiring a complete dataset based on a specified number of variables.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Raudenbush (2009). Analyzing effect sizes: Random-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 295-316). New York: Russell Sage Foundation.  
}
\examples{
# Sample data from Data Set III (Raudenbush, 2009; pp. 300)

require(MAc)
data(Data_ch16.1)

# Example

ComplData(temp, temp$mod1, temp$mod2, predictors= 2)
}
\seealso{
\code{\link{MAreg2}}
}
\keyword{data}

