\name{ComplData}
\alias{ComplData}
\title{Complete Dataset (by Moderator(s))
}
\description{ Function to reduce data set with complete data for x number of predictors in a multivariate moderator analysis. Will accept a data.frame with multiple rows per study and reduce the data.frame to one effect size per study based on recommended aggregation procedures for correlations (see chapter 10, Hunter & Schmidt, 2004; pp. 435-8).
}
\usage{
ComplData(meta, mod1, mod2 = NULL, mod3 = NULL, mod4 = NULL,
 mod5 = NULL, predictors = 1)
}
\arguments{
  \item{meta}{ data.frame with r (correlation coefficients) and n (sample size) for each study.
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
  \item{predictors}{ Number of moderators specified in function (to keep in reduced data.frame).
}
}
\value{ Reduced data.frame (with complete data) for each moderator entered into the function. This is primarily used as a convenience function for conducting multivariate meta-regressions.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\examples{
# Sample data from Data Set III (Cooper et al., 2009; pp. 265 & 587)

require(MAc)
data(Data_ch16.1)

# Example

ComplData(temp, temp$mod1, temp$mod2, predictors= 2)
}
\seealso{
\code{\link{MAreg2}}
}
\keyword{data}

