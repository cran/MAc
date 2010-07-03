\name{plotcon}
\alias{plotcon}
\title{Meta Regression Scatterplot
}
\description{Outputs a scatterplot from a fixed or random effects meta regression (continuous and/or categorical). 
}
\usage{
plotcon(es, var, mod, data, method= "random", modname=NULL, 
  title=NULL, ylim=c(0, 1), ...)
}
\arguments{
 \item{es}{r or z' effect size.
}
\item{var}{Vaiance of g.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
} 
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}. 
}
\item{data}{\code{data.frame} with values above.
}
  \item{modname}{Name of moderator to appear on x axis of plot. Default is NULL.
}
  \item{title}{Plot title. Default is NULL.
}
  \item{ylim}{Limits of y-axis with the first argument for the minimum y-value and the second for the maximum y value. Default is \code{c(0, 1)}.
}
\item{...}{ Additional arguments to be passed to ggplot.
  } 
}
\value{ Scatterplot with fixed or random effects regression line with size of visual points based on study weights, where the more precise studies have larger points. The ggplot2 package outputs the rich graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta analysis} (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{mareg}},
\code{\link{plotcat}}
}
\examples{
id<-c(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
.51,.40,.34,.42,.16)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
dat<-data.frame(id,n,r,mod1)
dat$var.r <- var_r(dat$r, dat$n) # MAc function to derive variance
dat$z <- r_to_z(dat$r)  # MAc function to convert to Fisher's z (z')
dat$var.z <- var_z(dat$n)  # MAc function for variance of z'
dat$mods2 <- factor(rep(1:2, 10))

# Example

\dontrun{plotcon(es = r, var = var.r, mod = mod1, data = dat, method= "random", 
modname= "Moderator") }
}
\keyword{ aplot }

