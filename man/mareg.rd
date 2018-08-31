\name{mareg}
\alias{mareg}
\title{ Meta-Regression 
}
\description{ Meta-regression function for a single or multiple predictor model. This function is a wrapper for the \code{rma()} function in the metafor package (Viechtbauer, W, 2010). Please see https://CRAN.R-project.org/package=metafor for details or for more advanced functionality with the \code{rma()} function. 
}
\usage{
mareg(formula, var, data, ztor = FALSE, method = "REML",  subset,  ...)
}
\arguments{
  \item{formula}{ This is a formula based function, similar to other functions in R (e.g., lm), where the criterion variable (e.g., r or z') is dependent on ('~') the predictor variables (e.g., moderators). The formula for two moderators would take this form: mareg(r ~ mod1 + mod2, var.r, data), where r is the criterion variable predicted by mod1 and mod2. The variance (var) of each r is var.r in this case.   
}
  \item{var}{ Variance of r or z'.
  }
\item{data}{Aggregated \code{data.frame} (see \code{agg} function for setting up the dataset for these analyses) with id, es (r or z'), var (variance of r or z') for each study. 
}
  \item{method}{ Default is \code{REML} (Restricted-Maximal Likelihood), which is the standard random effects method. For fixed effects, use \code{FE}. Other options are specified in the \code{metafor} package manual ('rma' function).
}
\item{ztor}{Default is FALSE. If TRUE, this assumes z' (Fisher's z) was used in the \code{es} argument and the analyist would like z' to be converted to r (for interpretive purposes) after analyzing in z'.
} 
 \item{subset}{ an optional vector specifying a subset of observations to be used in the fitting process.
  }
 \item{...}{ Additional arguments to be passed to rma().
  }  
}
\details{See Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. for the details  of the \code{rma()}  function. https://CRAN.R-project.org/package=metafor
}
\value{
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
\item{Pr(>|z|)}{ p-value (significance level).
}
\item{QE}{ Q-error. Measure of error in the model.
}
\item{QE.df}{ Degrees of freedom for Q-error.
}
\item{QEp}{ Q-error p-value (for homogeneity).
}
\item{QM}{ Q-model. Measure of model fit.
}
\item{QM.df}{ Degrees of freedom for Q-model.
}
\item{QMp}{ Q-between p-value (for homogeneity). QM and QMp provide the test of whether the moderator variable(s) account for significant variance among effect sizes.
} 
}
\references{Wolfgang Viechtbauer (2010). metafor: Meta-Analysis Package for
  R. R package version 1.1-0. https://CRAN.R-project.org/package=metafor
}
\seealso{
\code{\link{wd}},
\code{\link{plotcon}}
}
\examples{
# install metafor
# install.packages('metafor', dependencies = TRUE)

# Sample data
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



# Examples

# Random Effects
mareg(r~ mod1 + mods2, var = var.r, method = "REML",  data = dat)

# Fixed Effects
mareg(r~ mod1 + mods2, var = var.r, method = "FE",  data = dat)  
}



