\name{MAc-package}
\alias{MAc-package}
\alias{MAc}
\docType{package}
\title{Meta-Analysis with Correlations
}
\description{ This package contains a variety of functions relevant to conducting a correlational meta-analysis using recommended procedures as described in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, and Valentine, 2009). The goal in creating this package was to provide user-friendly functions to assist researchers in the process of conducting a meta-analysis, from the initial to final stages of their analytic endeavor. The meta-analyst can begin their project by using \code{MAc} functions to derive correlation coefficents from the primary studies, when statistics other than correlations are reported (e.g., t-test, p-value, or Cohen's d). Then, the analyst can aggregate all within-study effect sizes (while accounting for within-study correlations among outcome measures and eliminating any dependencies in the dataset) based on recommended procedures by Hunter & Schmidt (2004), calculate omnibus effect sizes under a fixed and random effects model, and assess for significant moderators (categorical and continuous, single and multi-predictor models) in the dataset. Finally, the meta-analyst can use one of several user-friendly graphics functions to visually represent their data. 
}
\details{
\tabular{ll}{
Package: \tab MAc\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2010-07-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The \code{MAc} package has integrated functions to facilitate the meta-analytic process at nearly every analytical stage. There are five broad areas of analysis that the \code{MAc} package targets:
 
1. Computations to Calculate Correlations: 

There are a variety of functions to compute r (correlation coefficients) from various designs reported in the primary studies. Most functions were derived from Borenstein's chapter (pp. 231-234) in \emph{The Handbook of Research Synthesis and Meta-Analysis} (Cooper, Hedges, & Valentine, 2009).
 For additional conversion formulas see the \bold{compute.es} package:  \url{https://CRAN.R-project.org/package=compute.es}


2. Within-Study Aggregation of Effect Sizes:

This package contains functions that have automated (i.e., will compute for all studies simultaneously) the process of aggregating within-study effect sizes while taking into account the dependencies among the within-study effect sizes (Hunter & Schmidt, 2004; Hedges & Olkin, 1985; Rosenthal et al., 2006). These functions fix the correlation between within-study effect sizes at .50 (Wampold et al., 1997) and will compute the correct aggregated effect size for all studies. \code{MAc} functions implement Hunter and Schmidt's (2004) recommendations for aggregating dependent correlations (see chapter 10, pp. 435-8).To our knowledge, this is the first statistical package/program that has explicitly utilized and automated this aggregation procedure, which has a dual effect of saving the researcher \bold{substantial} time while increasing the accuracy of their analyses. 

3. Fixed and Random Effects Omnibus Analysis: 

This package contains all the relevant functions to calculate fixed and random effects omnibus effect sizes, outputting the omnibus (i.e., overall) effect size, variance, standard error, upper and lower confidence intervals, and the Q-statistic (heterogeneity test). 

4. Moderator Analyses:

There are user-friendly functions to compute fixed and random effects moderator analyses. These include single and multiple predictor models for both categorical and continuous moderator data. 

5. Graphics:

This package has a variety of functions visually representing data. This includes boxplots and meta-regression scatterplots.

6. Sample of Additional Functions:

Export MA output to nicely formatted Word tables.

}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). \emph{The handbook of research synthesis and meta-analysis} (2nd edition). New York: Russell Sage Foundation.

Hunter, J. E., Schmidt, F. L. (2004). \emph{Methods of meta-analysis} (2nd edition). Thousand Oaks, CA: Sage. 

Wampold, B. E., Mondin, G. W., Moody, M., Stich, F., Benson, K., & Ahn, H. (1997). A meta-analysis of outcome studies comparing bona fide psychotherapies: Empiricially, 'all must have prizes.' \emph{Psychological Bulletin, 122(3)}, 203-215.
}
\keyword{ package }
\examples{ 
# Examples for each broad area:

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

# Note: for all the examples in this manual, we have made up data and manually 
# created variables and datasets. If conducting your own meta-analysis, 
# a more convenient way for using the functions is to import your data from
# a .csv file (with relevant variables in the dataset). One way to do this:
# dat <- read.csv(file.choose(), header = TRUE)
#
# Then, you can run the functions with this dataset and you do not need to 
# manually create your dataset, as we have done above.
  
  

# 1. Computations to Calculate Correlations:
 
# For example, suppose the primary study reported a t-test value for differences 
# between 2 groups and the total sample size. Then, running:

r_from_t (1.74, 30)  

# reported t-value (1.74) and sample size (30) will output the 
# correlation desired for use in the meta-analysis.

# 2. Within-Study Aggregation of Effect Sizes: 

agg(id = id, r = r, n = n, data=dat) 

# where data = data.frame with columns for id, r (correlation coefficient),
# and n (sample size) with multiple rows within each study (multiple
# correlations reported for each study). Outputs an aggregated data.frame 
# with 1 effect size per study. 

# 3. Fixed and Random Effects Omnibus Analysis

omni(es = z, var = var.z, data = dat, type="weighted", method = "random", ztor = TRUE)

# where data = data.frame with columns for id, es (r or z')
# , var (variance of r or z'), n (sample size). ztor = if using z', should
# it be converted back to r? see omni documentation for more details.
 
# 4. Moderator Analyses:

# Random effects
mareg(z~ mod1 + mods2, var = var.z, method = "REML", ztor = TRUE, data = dat) 
 
# where data = data.frame with columns for es (r or z'),
# var (variance of r or z') and moderators.

# 5. Graphics:

\dontrun{plotcon(es = r, var = var.r, mod = mod1, data = dat, method= "random", 
modname= "Moderator") }

# Additional Functions

# Export MA output to nicely formatted Word tables.

# install R2wd
# install.packages('R2wd', dependencies = TRUE)

# mareg fuction

temp <- mareg(z~ mod1 + mods2, var = var.z, method = "REML", ztor = TRUE, data = dat) 

# Export data to Word in formatted table

#  wd(temp, get = TRUE, new = TRUE)
}
\seealso{
\bold{RcmdrPlugin.MAc} package:  \url{https://CRAN.R-project.org/package=RcmdrPlugin.MAd};

\bold{MAd} package:  \url{https://CRAN.R-project.org/package=MAc};

\bold{RcmdrPlugin.MAd} package:  \url{https://CRAN.R-project.org/package=RcmdrPlugin.MAc};

\bold{compute.es} package:  \url{https://CRAN.R-project.org/package=compute.es};

\bold{metafor} package:  \url{https://CRAN.R-project.org/package=metafor}
}
