\name{plotcat}
\alias{plotcat}
\title{Categorical Moderator Graph
}
\description{ Outputs a rich and detailed boxplot graphic for each level of the specified moderator (under a fixed or random effects model). 
}
\usage{
plotcat(es, var, mod, data,  method="random",  modname=NULL,  title=NULL, ...)
}
\arguments{
 \item{es}{r or z' effect size.
}
\item{var}{Vaiance of es.
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
\item{...}{ Additional arguments to be passed to ggplot.
  } 
}
\value{ Boxplot graph with median, max, min, and outliers from a fixed or random effects categorical moderator analysis. Places jitter points (for each study) on the boxplots. The size of each point (representing a study in the analysis) are based on study weights where more precise studies have larger points. The ggplot2 package outputs the graphics.
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{macat}},
\code{\link{plotcon}}
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

\dontrun{plotcat(es = r, var = var.r, mod = mods2, data = dat, method= "random",
 modname= "Moderator") }
}   
\keyword{aplot}

