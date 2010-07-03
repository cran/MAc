\name{macatC}
\alias{macatC}
\title{Direct Categorical Moderator Comparison 
}
\description{Function for a planned comparison between two levels of a moderator under a fixed or random effects model.
}
\usage{
macatC(x1, x2, es, var, mod, data, method= "random", type= "post.hoc", ztor = FALSE)
}
\arguments{
 \item{x1}{One level of categorical moderator.
}
  \item{x2}{Comparison level of same categorical moderator.
}
 \item{es}{r or z' effect size.
}
\item{var}{Variance of es.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
} 
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}.
 }
  \item{type}{\code{post.hoc} assumes the comparison was not planned prior to conducting the meta analysis. The a priori option, \code{planned}, assumes the researcher planned to conduct the analysis a priori. Default is \code{post.hoc} using the Scheffe post hoc statistical method.
  } 
\item{ztor}{Default is FALSE. If TRUE, this assumes z' (Fisher's z) was used in the \code{es} argument and the analyist would like z' to be converted to r (for interpretive purposes) after analyzing in z'.
} 
 \item{data}{\code{data.frame} with values above.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{

\item{diff}{ Mean difference between the two levels.
} 
\item{var.diff}{ Variance of diff.
}
\item{p}{ Significance level.
}
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{macat}},
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
dat$mods2 <- factor(rep(1:4, 5))
dat

# Example
macatC(1, 2, es=r, var=var.r, mod=mods2, data=dat,  method= "random", 
type= "post.hoc", ztor = FALSE) 


}
\keyword{models}
