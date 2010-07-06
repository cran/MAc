\name{macat}
\alias{macat}
\title{Categorical Moderator Analysis 
}
\description{Computes single predictor categorical moderator analysis under a fixed or random effects model.
}
\usage{
macat(es, var, mod, data, method= "random", ztor = FALSE)
}
\arguments{
 \item{es}{r or z' effect size.
}
\item{var}{Variance of es.
}
  \item{mod}{Categorical moderator variable used for moderator analysis.
} 
 \item{method}{ Default is \code{random}. For fixed effects, use \code{fixed}. 
}
 \item{ztor}{Default is FALSE. If TRUE, this assumes z' (Fisher's z) was used in the \code{es} argument and the analyist would like z' to be converted to r (for interpretive purposes) after analyzing in z'.
}
 \item{data}{\code{data.frame} with values above.
}
}
\details{See Konstantopoulos & Hedges (2009; pp. 280-288) for the computations used in this function.
}
\value{

\item{mod}{ Level of the categorical moderator.
} 
\item{k}{ Number of studies for each level of the moderator.
}
\item{estimate}{ Mean effect size of each level of the moderator.
}
\item{ci.l}{ Lower 95\% confidence interval.
}
\item{ci.u}{ Upper 95\% confidence interval.
}
\item{z}{ z-score (standardized value).
}
\item{p}{ Significance level.
}
\item{var}{ Variance of effect size.
}
\item{se}{ Square root of variance.
}
\item{Q}{ Q-statistic (measure of homogeneity).
}
\item{df}{ Degrees of freedom for Q-statistic.
}
\item{p.h}{p-value for homogeneity within that level of the moderator. 
}
\item{I2}{ Proportion of total variation in effect size that is due to heterogeneity rather than chance (see Shadish & Haddock, 2009; pp. 263).
}
\item{Q}{ Q-statistic overall. Note: Whether fixed or random effects analyses are conducted, the Q statistic reported is for the fixed effect model. Therefore, Qb + Qw != Q in the random effects output.
} 
\item{Qw}{ Q-within (or error). Measure of within-group heterogeneity.
}
\item{Qw.df}{ Degrees of freedom for Q-within.
}
\item{Qw.p}{ Q-within p-value (for homogeneity).
}
\item{Qb}{ Q-between (or model). Measure of model fit.
}
\item{Qb.df}{ Degrees of freedom for Q-between.
}
\item{Qb.p}{ Q-between p-value (for homogeneity). Qb and Qb.p provide the test of whether the moderator variable(s) account for significant variance among effect sizes.
} 
}
\references{Konstantopoulos & Hedges (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 279-293). New York: Russell Sage Foundation.  

Shadish & Haddock (2009). Analyzing effect sizes: Fixed-effects models. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta analysis} (pp. 257-278). New York: Russell Sage Foundation. 
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{plotcat}},
\code{\link{wd}}
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
dat

# Example

# Random effects
macat(es = z, var= var.z, mod = mods2, data = dat, ztor = TRUE, method= "random")
}
\keyword{models}
