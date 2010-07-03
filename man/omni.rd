\name{omni}
\alias{omni}
\title{Omnibus Effect Size (Fixed and Random Effects) 
}
\description{Computes fixed and random effects omnibus effect size for correlations. 
}
\usage{
omni(es, var, data, type="weighted", method = "random", ztor = FALSE)
}
\arguments{
\item{es}{r or z' effect size.
}
\item{var}{Variance of es.
}
 
  \item{type}{\code{weighted} or \code{unweighted}. Default is \code{weighted}. Use the \code{unweighted} variance method only if Q is rejected and is very large relative to the number of studies in the meta-analysis. 
}
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}. 
}
\item{ztor}{Default is FALSE. If TRUE, this assumes z' (Fisher's z) was used in the \code{es} argument and the analyist would like z' to be converted to r (for interpretive purposes) after analyzing in z'.
} 
\item{data}{\code{data.frame} with above values.
}
}
\value{
Fixed and random effects:
 
\item{k}{ Number of studies in the meta-analysis.
}
\item{estimate}{ Unstandardized regression coefficient estimate.
} 
\item{se}{ Standard error of the estimate coefficient.
}
\item{z}{ z-value.
}
\item{ci.l}{ Lower 95\% confidence interval.
}
\item{ci.u}{ Upper 95\% confidence interval.
}
\item{p}{ Significance level.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df.Q}{ Degrees of freedom for Q-statistic.
}
\item{Qp}{ Q-statistic p-value (assesses overall homogeneity between studies).
}
\item{I2}{ Proportion of total variation in effect size that is due to systematic differences between effect sizes rather than by chance (see Shadish & Haddock, 2009; pp. 263).
}
}
\references{ Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
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

omni(es = z, var = var.z, data = dat, type="weighted", method = "random", ztor = TRUE)
}
\keyword{ models }

