\name{ComplData}
\alias{ComplData}
\title{Complete Dataset)
}
\description{ Function to reduce data set with complete data for 1 predictor. Will accept a \code{data.frame} with multiple rows per study and reduce the \code{data.frame} to one effect size per study based on recommended aggregation procedures for correlations (Hunter & Schmidt, 2004; pp. 435-8).
}
\usage{
ComplData(meta, mod)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study. The function requires that the column name for the correlation coefficient is 'r' and the sample size is 'n'.
}
  \item{mod}{ Moderator variable wanting to be kept for further analysis. 
}
}
\value{ Reduced \code{data.frame} (with complete data) for the moderator entered into the function. 
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

ComplData(temp, temp$mod1)
}
\keyword{data}

