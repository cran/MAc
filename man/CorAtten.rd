\name{CorAtten}
\alias{CorAtten}
\title{Correction for Attenuation
}
\description{ Used to correct for attenuated effect sizes due to measurement unreliability.
}
\usage{
CorAtten(meta, xx, yy)
}
\arguments{
  \item{meta}{ \code{data.frame} with r (correlation coefficients) and n (sample size) for each study.
}
  \item{xx}{ Column for reliability of predictor variable ("independent variable"). 
}
  \item{yy}{ Column for reliability of outcome variable ("dependent variable").
}
}
\value{ \code{data.frame} with a new column for correlations corrected for measurement unreliability (\code{r.adj}), updated standard errors, variance, and study weights based on these corrected values (see Hunter & Schmidt, 2004; pp. 97-98). Studys without reliability information will remain unchanged, as will their standard errors, variances and weights.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{Hunter, J. E., Schmidt, F. L. (2004). \emph{Methods of meta-analysis} (2nd edition). Thousand Oaks, CA: Sage. 
}
\examples{
# Sample data:

id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
xx<-c(.88,.86,.83,.64,.89,.84,.89,.83,.99,.88,.81,.77,.86,.70,.79,
 .71,.80,.74,.82,.86)  # Reliability of "independent variable"
yy<-c(.99,.86,.83,.94,.89,.94,.89,.93,.99,.88,.81,.77,.86,.70,.79,
 .71,.80,.94,.92,.96)  # Reliability of "dependent variable"
   
data<-data.frame(id,n,r,xx,yy)

# Example        
CorAtten(data,data$xx,data$yy) 
}
\keyword{data}

