\name{MAc-package}
\alias{MAc-package}
\alias{MAc}
\docType{package}
\title{Meta-Analysis for Correlations
}
\description{ This package contains all the relevant functions to conduct a correlational meta-analysis using standard procedures as described in Cooper, Hedges, and Valentines Handbook of Research Synthesis and Meta-Analysis (2009). The goal in creating this package was to assist researchers from the initial to final stages of conducting a research synthesis using the most up to date statistical procedures. The meta-analyst can begin their project by using functions in this package to transform reported statistics from primary studies into correlation coefficients.  
}
\details{
\tabular{ll}{
Package: \tab MAc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-02-02\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
This package has integrated functions to facilitate the meta-analytic process (for correlations) at every analytical stage. There are five broad areas of analysis that the MAc package targets:
 
1. Computations to Calculate Correlations: 

There are a variety of functions to compute r (correlation coefficients) from various designs reported in the primary studies. Most functions were derived from section 12.4 & Table 12.4 (pp. 231-234) of Cooper, Hedges, and Valentines Handbook of research synthesis and meta-analysis (2009).

2. Within-Study Aggregation of Effect Sizes:

Automated functions (i.e., will compute for all studies simultaneously) for aggregating within-study effect sizes (which account within study dependencies of effect sizes). These functions fix the correlation between predictor variables at .50 and will compute the correct aggregated effect size for all studies. Function implements Hunter and Schmidts (2004) approach to aggregation of dependent correlations (see chapter 10, pp. 435-8).

3. Fixed and Random Effects Omnibus Analysis: 

This package will calculate fixed and random effects omnibus effect sizes for correlations, outputting the omnibus (i.e., overall) effect size, variance, standard error, upper and lower confidence intervals, and the Q-statistic (heterogeneity test). 

4. Moderator Analyses:

Functions to compute fixed and random effect moderator analyses. This includes single and multiple predictor models for both categorical and continuous data.

5. Graphics:

Functions for a variety of meta-analytic graphics, including standard forest plots, funnel plots, boxplots, meta-regression scatterplots, and unique graphics, such as multi-predictor meta-regression graphs.
}
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\references{ Cooper, H., Hedges, L.V., & Valentine, J.C. (2009). The handbook of research synthesis and meta-analysis (2nd edition). New York: Russell Sage Foundation.

Hunter, J. E., Schmidt, F. L. (2004). Methods of meta-analysis (2nd edition). Thousand Oaks, CA: Sage. 
}
\keyword{ package }
\examples{ 
# Examples for each broad area:

id<-rep(1:20)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
 .51,.40,.34,.42,.16)
mod1<-c(1,2,3,4,1,2,8,7,5,3,9,7,5,4,3,2,3,5,7,1)
data<-data.frame(id,n,r,mod1)

#Computations to Calculate Correlations:
 
# Primary study reported t (t-test value of differences between 2 
# groups), n (total sample size)

r_from_t (1.74, 30)  

# reported t-value (1.74) and sample size (30) 

# Within-Study Aggregation of Effect Sizes: 

MetaR(data)  

# where data = data.frame with columns for r and n with multiple rows per # study. Outputs an aggregated data.frame with 1 effect size per study 

#Fixed and Random Effects Omnibus Analysis

OmnibusES(data) 

# where data = data.frame with columns for r and n
 
# Moderator Analyses:

MAreg1(data, data$mod1) 
 
# where data = data.frame with columns for r and n data$mod1 = moderator # variable of interest.

# Graphics:

MAregGraph(data, data$mod1)  

# same as above but will output a meta-regression graph
}
