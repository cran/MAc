##==================             MAc             ================##      
##==================  Meta-Analysis Correlations ================##

# updated: 03.09.10

# Package created by AC Del Re & William T. Hoyt
# This package contains all the relevant functions to conduct a
# correlational meta-analysis using standard procedures as
# described in Cooper, Hedges, & Valentine's Handbook of
# Research Synthesis and Meta-Analysis (2009).
# This package requires the 'ggplot2' package, 'plyr' package.

require('ggplot2')
require('plyr')
#require('psych') 

##=== New Functions ===##
# 2.24.10: added (1) conversion to d functions, (2) CatComp
# update
# 2.23.10: rm agg_r2 (problematic) & CatMod update
# new element of function for mods 2.20.10. Big thanks to 
# Jim Holtman for the function to randomly select mods
# integrated into agg_r2()
# New function: 2.16.10
# Multiple predictor regression model fit
# (required inputs are lm saved as objects)

MRfit <- function( ...) {
  models <- list(...)
  fit<- do.call(anova, models)
  fit$R2 <- 1 - (fit$RSS / fit$RSS[1])
  fit$R2.change <- c(NA, diff(fit$R2))
  fit$R2[1] <- NA
  names(fit) <- c("df.Q", "Qw", "predictors", "Qb", "F",  "Pvalue", 
                  "R^2", "R^2.change")
  return(fit)
}



# a formula for attenuation already created:  correct.cor(x, y)


##=== Preliminary Steps ===##

# Import data into R:
# 1. Save main data file (excel or spss) to .csv [e.g.,  see save options in excel]
# 2. Import .csv file into R by setting the working directory to the location of 
#    your data file,  e.g.:
#    setwd("C:/Users/User/Documents/TA Meta-Analy/Horvath_2009/ANALYSIS/12-10-09") 
#    and then import data,  e.g.:
#    data <- read.csv("Alliance_1-30-10.csv", header=TRUE, na.strings="") 

##==== Data Manipulation ====##


# set numeric variables to numeric, e.g.:
# data$r <- as.numeric(as.character(data$r))

# set categorical variables to factors or character,  e.g.:
# data$id <- as.character(data$id)

# fix data with errors in factor names, requires car package,e.g.:
# library(car)
# data$outcome3 <- recode(data$outcome2, 'c("?",  "adherence", "compliance", "depression", 
#                         "depression ", "wellbeing", "work", "GAS")="Other"; 
#                         c("GSI", "SCL", "BSI")="SCL"; c("dropout")="Dropout"; 
#                         else= "Other"')    

##============ COMPUTATIONS TO CALCULATE EFFECT SIZES ================##

# Formulas for computing d, var(d), g (bias removed from d), and var(g)
# in designs with independent groups.
# Section 12.3.1 & Table 12.1 (Cooper et al., 2009; pp. 226-228)

# Computing d and g, independent groups
# (12.3.1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# sd.1 (treatment standard deviation at post-test), sd.2 (comparison 
# standard deviation at post-test), n.1 (treatment), n.2 (comparison/control).

d_to_g <- function(d, var.d, n.1, n.2) {
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  out<-cbind(g, var.g)
  return(out)
}

mean_to_d <- function(m.1,m.2,sd.1,sd.2,n.1, n.2) {
  s.within<-sqrt(((n.1-1)*sd.1^2+(n.2-1)*sd.2^2)/(n.1+n.2-2))
  d<-(m.1-m.2)/s.within
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var.d)
  return(out)
}

# (1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control).

mean_to_d2 <- function(m.1,m.2,s.pooled,n.1, n.2) {
  d<-(m.1-m.2)/s.pooled
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (2) Study reported: 
# t (t-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).

t_to_d <- function(t, n.1, n.2) {
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (3) Study reported: 
# f (F-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).


f_to_d <- function(f,n.1, n.2) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (4) Study reported: 
# p-value (for ONE-tailed test), n.1 (treatment), n.2 (comparison/control).

p_to_d1 <- function(p, n.1, n.2) {
  pxtwo<-p*2
  df<-(n.1+n.2)-2  
  TINV<-qt((1-pxtwo/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (5) Study reported: 
# p-value (for TWO-tailed test), n.1 (treatment), n.2 (comparison/control).

p_to_d2 <- function(p, n.1, n.2) {
  df<-(n.1+n.2)-2  
  TINV<-qt((1-p/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
  var_d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# rtod converts Pearson r to Cohen's d (also retains (N-1)/N).
r_to_d <- function(r, N) { 2*r*sqrt((N-1)/(N*(1-r^2)))*abs(r)/r    }

# Formulas for computing d and var(d) in designs with independent groups
# using ANCOVA. Section 12.3.3 & Table 12.3 (Cooper et al., 2009; pp. 228-230).

# Computing d and g from ANCOVA
# (12.3.3) Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# sd.adj (adjusted standard deviation), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

ancova_to_d1 <- function(m.1.adj,m.2.adj,sd.adj,n.1, n.2, R, q) {
  s.within<-sd.adj/sqrt(1-R^2)
  d<-(m.1.adj-m.2.adj)/s.within
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}


# Table 12.3 (1) Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

ancova_to_d2 <- function(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q) {
  d<-(m.1.adj-m.2.adj)/s.pooled
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (2) Study reported: 
# t (t-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

tt.ancova_to_d <- function(t, n.1, n.2, R, q) {
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (3) Study reported: 
# f (F-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control),R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

f.ancova_to_d<-function(f,n.1, n.2, R, q) {
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (4) Study reported: 
# p-value (for ONE-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

p.ancova_to_d1 <- function(p, n.1, n.2, R, q) {
  pxtwo<-p*2
  df<-(n.1+n.2)-2  
  TINV<-qt((1-pxtwo/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# (5) Study reported: 
# p-value (for TWO-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

p.ancova_to_d2 <- function(p, n.1, n.2, R, q) {
  df<-(n.1+n.2)-2  
  TINV<-qt((1-p/2),df)
  d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var_d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  out<-cbind(d,var_d)
  return(out)
}

# computing d from odds ratio

or_to_d <- function(or) {
  lor <- log(or)
  d <- lor*sqrt(3)/pi
  out <- cbind(lor, d)
  return(out)
}  

# computing d from log odds ratio

lor_to_d <- function(lor, var.lor) {
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  out <- cbind(d, var.d)
  return(out)
}  

# compute or from proportions

prop_to_or <- function(p1, p2, n.ab, n.cd) {
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  out <- cbind(or,lor,var.lor)
  return(out)
}

prop_to_d <-function(p1, p2, n.ab, n.cd) {
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  out <- cbind(or,lor,var.lor, d, var.d)
  return(out)
}

# Odds Ratio to d: if have info for 'failure' in both conditions 
# (B = # tmt failure; D = # non-tmt failure) and the sample size
# for each group (n.1 & n.2 respectively):

fail_to_d <- function(B, D, n.1, n.0) {
  A <- n.1 - B  # tmt success
  B <- B        # tmt failure
  C <- n.0 - D  # non-tmt success
  D <- D        # non-tmt failure
  p1 <- A/n.1   # proportion 1 
  p2 <- C/n.0   # proportion 2
  n.ab <-  A+B  # n of A+B
  n.cd <-  C+D  # n of C+D        
  or <- (p1 * (1 - p2))/(p2 * (1 - p1))  # odds ratio
  lor <- log(or)  # log odds ratio
  var.lor <-  1/A + 1/B + 1/C + 1/D  # variance of log odds ratio
  #var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor * sqrt(3)/pi  # conversion to d
  var.d <- 3 * var.lor/pi^2  # variance of d
  out <- cbind(or, lor, var.lor, d, var.d)
  return(out)
}

##============ COMPUTATIONS TO CALCULATE CORRELATIONS ================##

# Formulas for computing r in designs with independent groups. 
# Section 12.4 & Table 12.4 (Cooper et al., 2009; pp. 231-234).

# (1) Computing variance of r,  z', and variance of z':

var_r <-  function(r, n) ((1 - r^2)^2)/(n - 1)  #calulate variance of r     
r_to_z  <-  function(r) 0.5*log((1 + r)/(1 - r))  #convert to z' 
var_z <- function(n) 1 / (n - 3)  #variance of z'

# Formulas for computing r in designs with independent groups. 
# Section 12.4 & Table 12.4 (Cooper et al., 2009; pp. 231-234).

# (1) Study reported: 
# t (t-test value of differences between 2 groups), n (total sample size)

r_from_t <- function(t, n) {
  r <- sqrt((t^2)/(t^2 + n-2))
  var_r <- ((1-r^2)^2)/(n-1)
  out <- cbind(r, var_r)
  return(out)
}

# Converting d (mean difference) to r where n.tmt = n.comparison 
# (Section 12.5.4; pp. 234)

r_from_d <- function(d,  var.d,  a=4) {
  r <- d/sqrt((d^2) + a)
  var_r <- (a^2*var.d)/(d^2 + a)^3
  out <- cbind(r, var_r)
  return(out)
}

# Converting d to r where n.tmt (not) = n.comparison (Section 12.5.4; pp. 234)

r_from_d1 <- function(d,  n.1, n.2,  var.d) {
  a <- ((n.1 + n.2)^2)/(n.1*n.2)
  r <- d/sqrt((d^2) + a)
  var_r <- (a^2*var.d)/(d^2 + a)^3
  out <- cbind(r, var_r)
  return(out)
  }

# Converting Chi-squared statistic with 1 df to r

r_from_chi <- function(chi.sq,  n) sqrt(chi.sq/n)

##============ WITHIN STUDY AGGREGATION OF EFFECT SIZES =============##

# Functions for aggregating within-study effect sizes (accounting for dependencies).
# Required inputs are a data.frame containing id (study id number) and 
# r (correlations).  This function fixes the correlation between predictor 
# variables at .50 and will compute the correct aggregated effect size for
# all studies. Function implements Hunter & Schmidt (2004)
# approach to aggregation of dependent r's (see chapter 10,  pp. 435-8).

aggrs <- function(r) {
  # Intermediate level function to compute weighted aggregation of effect sizes 
  k <- length(r)
   rbar <- .50
   r.xY <- sum(r)/(1 * sqrt(k  +  k*(k - 1) * rbar))
   return(r.xY)
}

# Intermediate level function to agg data and create new data.frame. 
# Required inputs are data.frame with r (correlation coefficients) and n (sample size)

agg_r <- function(meta) {
  agg <- with(meta,  aggregate(list(r = r),  by = list(id = id), aggrs))
  ns <- with(meta,  aggregate(list(n = n),  by = list(id = id),  mean))
  ns$n <- round(ns$n, 0)
  agg.data  <-  merge(agg,  ns,  by='id') 
  return(agg.data)
}

# used to randomly select one level of mod if study reports multi levels

agg_r2 <- function(meta, mod) {
  agg <- with(meta,  aggregate(list(r = r),  by = list(id = id, mod=mod),  aggrs))
  ns <- with(meta,  aggregate(list(n = n),  by = list(id = id, mod=mod),  mean))
  meta  <-  merge(agg,  ns,  by= c('id','mod')) 
  meta <- as.data.frame(meta)
  meta$n <- round(ns$n, 0)
  return(meta)
  }

##=== Add Fixed and Random Effects Weights ===##
# Required input is a data.frame with column names id (study id), 
# r (correlation coefficient),  and n (sample size).
 
MetaR <-  function(meta) {  
  # Computes within study aggregation and adds fixed and random effects weights 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for 
  #   each study.
  # Returns:
  #   Aggregated effect size (one row per study) based on Hunter and Schmidt's (2004)
  #   approach to aggregation of dependent r's (see chapter 10,  pp. 435-8).
  #   Adds study weights based of a Fisher's z transformation of r's (see chapter 14, 
  #   Cooper et al., 2009)
  meta <- agg_r(meta)
  meta$var.r <- var_r(meta$r, meta$n)
  meta$var.r <-  ifelse(is.na(meta$var.r), ((1-meta$r^2)^2)/(meta$n-1), meta$var.r)
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))  #computing r to z' for each study
  meta$var.z <- 1/(meta$n-3)  # computing var.z for each study
  meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
  meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
  meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
  meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)  
  meta$z.score <- meta$z/sqrt(meta$var.z)
  meta$p.value <- 2*(1-pt(abs(meta$z.score),  meta$n-1))
  meta$wi <-  1/meta$var.z  # computing weight for each study
  meta$wiTi <- meta$wi*meta$z  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$z)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$r))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + meta$var.z  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  meta$l.ci95z <- NULL 
  meta$u.ci95z <- NULL
 return(meta)
}

##================= FIXED AND RANDOM EFFECTS OMNIBUS ===============##
# Function to calculate fixed and random effects omnibus effect size for correlations,  
# outputing omnibus effect size,  variance,  standard error,  upper and lower 
# confidence intervals,  and heterogeneity test.
# Required inputs are a data frame with id (study id) and r (correlation 
# coefficient)variable names.

OmnibusES<-  function(meta,  var="weighted" ) {
  # Computes fixed and random effects omnibus effect size for correlations.  
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   var:  "weighted" or "unweighted". "weighted" is the default. Use the 
  #   unweighted variance method only if Q is rejected and is very large relative to k.   
  # Returns:
  #   Fixed and random effects omnibus effect size, variance, standard error, 
  #   upper and lower confidence intervals, p-value, Q (heterogeneity test), I2
  #   (I-squared--proportion of total variation in tmt effects due to heterogeneity 
  #   rather than chance). 
  meta <- MetaR(meta)
  k <- length(!is.na(meta$r)) # number of studies
  df <- k-1 
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate omnibus
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate omnibus
  Tz.agg <- sum.wiTi/sum.wi  # omnibus z' 
  var.Tz.agg <- 1/sum.wi  # omnibus var.z
  se.Tz.agg <- sqrt(var.Tz.agg) 
  z.value  <-  Tz.agg/se.Tz.agg
  p.value <-  round(2*(1-pt(abs(z.value),  df)), 6)  
  T.agg <- (exp(2*Tz.agg)-1)/(exp(2*Tz.agg) + 1)  # fixed effect (FE) omnibus effect size (r)
  var.T.agg <- (exp(2*var.Tz.agg)-1)/(exp(2*var.Tz.agg) + 1)  # FE omnibus var.r
  lower.ci.Tz <- Tz.agg-1.96*se.Tz.agg
  upper.ci.Tz <- Tz.agg + 1.96*se.Tz.agg
  lower.ci <- (exp(2*lower.ci.Tz)-1)/(exp(2*lower.ci.Tz) + 1)  # FE lower CI
  upper.ci <- (exp(2*upper.ci.Tz)-1)/(exp(2*upper.ci.Tz) + 1)  # FE upper CI
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  # FE homogeneity test
  I2 <- (Q-(k-1))/Q  # I-squared 
  I2 <- ifelse(I2<0, 0, I2)        
  I2 <- paste(round(I2*100, 4),  "%",  sep="")                        
  p.homog <- pchisq(Q, df, lower=FALSE)  # <.05 = sig. heterogeneity  
  # random effects #
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
  sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
  sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
  Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
  if(var == "weighted") {
    var.Tz.agg.tau <-  1/sum.wi.tau 
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    z.valueR  <-  Tz.agg.tau/se.Tz.agg.tau
    p.valueR <-  round(2*(1-pt(abs(z.valueR),  df)), 6)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    lower.ci.Tz.tau <- Tz.agg.tau-1.96*se.Tz.agg.tau
    upper.ci.Tz.tau <- Tz.agg.tau + 1.96*se.Tz.agg.tau
    lower.ci.tau <- (exp(2*lower.ci.Tz.tau)-1)/(exp(2*lower.ci.Tz.tau) + 1)
    upper.ci.tau <- (exp(2*upper.ci.Tz.tau)-1)/(exp(2*upper.ci.Tz.tau) + 1)
  }
  if(var == "unweighted") {  # unweighted variance method
    var.agg <- (sum(meta$z^2)-sum(meta$z)^2/k)/(k-1) #14.20
    q.num <- (1/k)*sum(meta$var.z)                                   
    unwgtvar.Tz.agg.tau <- var.agg-q.num  #14.22
    var.Tz.agg.tau <-  ifelse(unwgtvar.Tz.agg.tau <= 0, 0, unwgtvar.Tz.agg.tau)  #if var < 0,  its set to 0
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    z.valueR  <-  Tz.agg.tau/se.Tz.agg.tau
    p.valueR <-  round(2*(1-pt(abs(z.valueR),  df)), 6)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    lower.ci.Tz.tau <- Tz.agg.tau-1.96*se.Tz.agg.tau
    upper.ci.Tz.tau <- Tz.agg.tau + 1.96*se.Tz.agg.tau
    lower.ci.tau <- (exp(2*lower.ci.Tz.tau)-1)/(exp(2*lower.ci.Tz.tau) + 1)  
    upper.ci.tau <- (exp(2*upper.ci.Tz.tau)-1)/(exp(2*upper.ci.Tz.tau) + 1)
  } 

  Fixed <- list(FixedEffects=c(k=k, r=T.agg,  var.r=var.T.agg,  se=se.Tz.agg, 
                l.ci=lower.ci,  u.ci=upper.ci,  z.value=z.value,  p.value=p.value,
                Q=Q, df.Q=df, p_homog=p.homog, I2=I2))
  Random <- list(RandomEffects=c(k=k, r=T.agg.tau, var.r=var.T.agg.tau,  
                 se=se.Tz.agg.tau,  l.ci=lower.ci.tau, u.ci=upper.ci.tau, 
                 z.value=z.valueR, p.value=p.valueR, Q=Q, df.Q=df,  p_homog=p.homog, 
                 I2=I2))
  omni.data <- as.data.frame(c(Fixed, Random))      
  omni.data$Omnibus <- c("k", "ES", "var.ES", "SE", "CI.lower", 
                         "CI.upper", "Z", "p", "Q", "df", "p.hetero", 
                         "I2")
  omni.data <- omni.data[c(3, 1, 2)]
  omni.data <- as.data.frame(omni.data)
  row.names(omni.data) <- NULL
  return(omni.data)
}

# Now,  if there is significant heterogeneity (p_homog < .05),  look for moderators.

##================= Categorical Moderator Analysis ================##
# requires plyr package

# Fixed effects intermediate level functions:

f1 <- function(meta) c(k=length(meta$k), sum.wi=sum(meta$wi, na.rm=TRUE), 
                       sum.wiTi=sum(meta$wiTi, na.rm=TRUE), sum.wiTi2=sum(meta$wiTi2,
                       na.rm=TRUE))
f2 <- function(modsum) c(K=modsum$k, ES=modsum$sum.wiTi/modsum$sum.wi, 
                         Var=1/modsum$sum.wi, 
                         Q=modsum$sum.wiTi2-(modsum$sum.wiTi^2)/modsum$sum.wi)

# pick_one added 2.17.10 (thanks to Hadley Wickham for assistance)

pick_one <- function(x) {
 if (length(x) == 1) return(x)
 sample(x, 1)
}

# Intermediate level function: sums relevant es information by each moderator level 
cat_sum1 <- function(meta, mod) {
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #          function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  mod <- meta$mod
  sum.mod <- by(meta, mod, f1)
  sums <- cbind(expand.grid(dimnames(sum.mod)), do.call(rbind, sum.mod))
  return(sums)
}
# Intermediate level function: adds an overall category to moderators

cat_sum2 <-  function(meta, mod) {
  catsum <- cat_sum1(meta, mod)
  wi <- sum(catsum$sum.wi)
  wiTi <- sum(catsum$sum.wiTi)
  wiTi2 <- sum(catsum$sum.wiTi2)
  k <- sum(catsum$k)
  mod <- catsum$mod
  catsums <- rbind(catsum, data.frame(mod="Overall", sum.wi=wi, sum.wiTi=wiTi, 
                   sum.wiTi2=wiTi2, k=k))
  return(catsums)
}

# Intermediate level function: derives mean ES & Var for each group

mod_mean <- function(meta,  mod) {
  modsum <- cat_sum2(meta,  mod)
  sum.values <- by(modsum, modsum$mod, f2)
  sums <- cbind(expand.grid(dimnames(sum.values)), do.call(rbind, sum.values))
  colnames(sums) <- c("Mod", "K", "ES", "Var", "Q")
  return(sums)
}

# Intermediate level function: adds 95% confidence intervals
mod_sig.z <- function(meta,  mod) {
  mod.es <- mod_mean(meta, mod)
  mod.es$L.95ci <- mod.es$ES-1.96*sqrt(mod.es$Var)
  mod.es$U.95ci <- mod.es$ES + 1.96*sqrt(mod.es$Var)
  return(mod.es)}

# Random effects intermediate level functions:

f3 <- function(meta) c(k=length(meta$k), sum.wi.tau=sum(meta$wi.tau, 
                       na.rm=TRUE),  sum.wiTi.tau=sum(meta$wiTi.tau, na.rm=TRUE), 
                       sum.wiTi2.tau=sum(meta$wiTi2.tau, na.rm=TRUE))
f4 <- function(modsum) c(K=modsum$k, ES=modsum$sum.wiTi.tau/modsum$sum.wi.tau, 
                         Var=1/modsum$sum.wi.tau,  
                         Q=modsum$sum.wiTi2.tau-(modsum$sum.wiTi.tau^2)/modsum$sum.wi.tau)

# Intermediate level function: sums relevant es information by each moderator level 

cat_sum1R <- function(meta, mod) {
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #          function(.data) .data[sample(nrow(.data), 1),]))
  meta$z <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <- 1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-(k - 1))/comp  #random effects variance
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  sum.mod <- by(meta, mod, f3)
  sums <- cbind(expand.grid(dimnames(sum.mod)), do.call(rbind, sum.mod))
  names(sums)[1] <- "mod"
  return(sums)
}

# Intermediate level function: adds an overall category to moderators

cat_sum2R <- function(meta, mod) {
  catsum <- cat_sum1R(meta, mod)
  wi.tau <- sum(catsum$sum.wi.tau)
  wiTi.tau <- sum(catsum$sum.wiTi.tau)
  wiTi2.tau <- sum(catsum$sum.wiTi2.tau)
  k <- sum(catsum$k)
  mod <- meta$mod
  catsums <- rbind(catsum, data.frame(mod="Overall", sum.wi.tau=wi.tau, 
  sum.wiTi.tau=wiTi.tau,  sum.wiTi2.tau=wiTi2.tau, k=k))
  return(catsums)
}


# Intermediate level function: derives mean ES & Var for each group

mod_meanR <- function(meta,  mod) {
  modsum <- cat_sum2R(meta,  mod)
  sum.values <- by(modsum, modsum$mod, f4)
  sums <- cbind(expand.grid(dimnames(sum.values)), do.call(rbind, sum.values))
  colnames(sums) <- c("Mod", "K", "ES", "Var", "Q")
  return(sums)
}

# Intermediate level function: adds 95% confidence intervals

mod_sig.zR <- function(meta,  mod) {
  mod.es <- mod_meanR(meta, mod)
  mod.es$L.95ci <- mod.es$ES-1.96*sqrt(mod.es$Var)
  mod.es$U.95ci <- mod.es$ES + 1.96*sqrt(mod.es$Var)
  return(mod.es)}

# Single predictor categorical moderator function (fixed effects)

CatModf <-  function(meta,  mod) {
  # Computes single predictor categorical moderator analysis. Computations derived from 
  # chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Fixed effects moderator means per group, k per group, 95% confidence intervals,
  #   z-value, p-value, variances, standard errors, Q, df(Q), and I-squared.
  mod.sig <- mod_sig.z(meta, mod)
  mod.sig$ES <- (exp(2*mod.sig$ES)-1)/(exp(2*mod.sig$ES) + 1)
  mod.sig$Var <- (exp(2*mod.sig$Var)-1)/(exp(2*mod.sig$Var) + 1)
  mod.sig$Std.Err <- sqrt(mod.sig$Var)
  k <- mod.sig[mod.sig$Mod=="Overall", "K"]  
  df <- k-1
  mod.sig$z.value <- mod.sig$ES/mod.sig$Std.Err
  mod.sig$p.value <- round(2*(1-pt(abs(mod.sig$z.value),  df)), 6)
  mod.sig$L.95ci <- (exp(2*mod.sig$L.95ci)-1)/(exp(2*mod.sig$L.95ci) + 1)
  mod.sig$U.95ci <- (exp(2*mod.sig$U.95ci)-1)/(exp(2*mod.sig$U.95ci) + 1)
  mod.sig$df.Q <- mod.sig$K-1
  mod.sig$p_homog <- ifelse(mod.sig$df.Q==0,  1,  round(pchisq(mod.sig$Q, mod.sig$df.Q, lower=FALSE), 5)) 
  mod.sig$I2 <- (mod.sig$Q-(mod.sig$K-1))/mod.sig$Q   #I-squared  
  mod.sig$I2 <- ifelse(mod.sig$I2<0, 0, mod.sig$I2)     
  mod.sig$I2 <- paste(round(mod.sig$I2*100, 4),  "%",  sep="")
  mod.sig2 <- mod.sig[, c(1, 2, 3, 6, 7, 9, 10, 4, 8, 5, 11, 12, 13)] #effects due to heterogeneity rather than chance (p263)
  rownames(mod.sig2) <- NULL
  colnames(mod.sig2) <- c("Mod", "k", "ES", "CI.lower", "CI.upper",  "Z", 
  "p",  "var.ES", "SE",  "Q", "df", "p.hetero",  "I2")
  return(mod.sig2)
}

# Fixed effect single predictor categorical moderator Q-statistic

CatModfQ <- function(meta,  mod) {
  # Computes fixed effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Fixed effects moderator Q-statistic, Q-within & between, df(Qw & Qb), and 
  #   homogeneity p-value within & between levels.
  mod.sig <- CatModf(meta, mod)
  k <- mod.sig$k[mod.sig$Mod=="Overall"]                           #number of studies
  levels <- length(mod.sig$Mod)-1
  Qb.df <- levels-1 
  Qw.df <- k-levels
  Q <- mod.sig$Q[mod.sig$Mod=="Overall"]         #overall heterogeneity Q
  Qw <- sum(mod.sig$Q[!mod.sig$Mod=="Overall"])  #overall within-group heterogeneity statistic
  Qw_p.value <- 1-pchisq(Qw, Qw.df)      
  Qb <- Q-Qw                                     #overall between-group heterogeneity
  Qb_p.value <- 1-pchisq(Qb, Qb.df)
  mod.Qstat <- data.frame(Q, Qw, Qw.df, Qw_p.value, Qb, Qb.df, Qb_p.value)
  names(mod.Qstat) <- c("Q", "Qw", "df.w",  "p.w", "Qb", "df.b", "p.b")
  return(mod.Qstat)
}
         
# Function for planned comparisons between 2 levels of moderator (fixed effects)

CatCompf <- function(meta,  mod, x1, x2,  method= "post.hoc1") {
  # Directly compares 2 levels of a categorical moderator using a fixed effects model. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc1". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  modsig <- mod_sig.z(meta, mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "Var"]
  x2.var <- modsig[modsig$Mod==com2, "Var"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta 
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$z~meta$mod,weights=meta$wi)
    fit<-TukeyHSD(fit)
  }
  if(method == "post.hoc2") {  # post-hoc comparison (Scheffe method)
    z <- g/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- g-1.96*sqrt(var)
    U.95ci <- g + 1.96*sqrt(var)
    fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
    names(fit) <- c("diff", "var.diff",  
                    "p", "CI.lower", "CI.lower" )
  }
  if (method == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- g^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- g-1.96*sqrt(var)
   U.95ci <- g + 1.96*sqrt(var)
   fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
                   "p", "CI.lower", "CI.lower" )
  }
  return(fit) 
}

# Function for planned comparisons between 2 levels of moderator (random effects)

CatModr <-  function(meta,  mod) {
  # Computes single predictor random effects categorical moderator analysis. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Random effects moderator means per group, k per group, 95% confidence intervals,
  #   z-value, p-value, variances, standard errors, Q, df(Q), and I-squared.
  mod.sig <- mod_sig.zR(meta, mod)
  mod.sig$ES <- (exp(2*mod.sig$ES)-1)/(exp(2*mod.sig$ES) + 1)
  mod.sig$Var <-(exp(2*mod.sig$Var)-1)/(exp(2*mod.sig$Var) + 1)
  mod.sig$Std.Err <- sqrt(mod.sig$Var)
  k <- sum(!is.na(mod)) 
  df <-  k-1
  mod.sig$z.value <- mod.sig$ES/mod.sig$Std.Err
  mod.sig$p.value <- round(2*(1-pt(abs(mod.sig$z.value),  df)), 6)
  mod.sig$L.95ci <- (exp(2*mod.sig$L.95ci)-1)/(exp(2*mod.sig$L.95ci) + 1)
  mod.sig$U.95ci <-  (exp(2*mod.sig$U.95ci)-1)/(exp(2*mod.sig$U.95ci) + 1)
  mod.sig$Q <- mod.sig$Q
  mod.sig$df.Q <- mod.sig$K-1
  mod.sig$p_homog <- ifelse(mod.sig$df.Q==0, 1, round(pchisq(mod.sig$Q, mod.sig$df.Q, lower=FALSE), 5)) 
  mod.sig$I2 <- (mod.sig$Q-(mod.sig$K-1))/mod.sig$Q   #I-squared  
  mod.sig$I2 <- ifelse(mod.sig$I2<0, 0, mod.sig$I2)     
  mod.sig$I2 <- paste(round(mod.sig$I2*100, 4),  "%",  sep="")
  mod.sig2 <- mod.sig[, c(1, 2, 3, 6, 7, 9, 10, 4, 8, 5, 11, 12, 13)]
  rownames(mod.sig2) <- NULL
  colnames(mod.sig2) <- c("Mod", "k", "ES", "CI.lower", "CI.upper",  "Z", 
  "p",  "var.ES", "SE",  "Q", "df", "p.hetero",  "I2")
  return(mod.sig2)
}  

# Intermediate function for random effects Q-statistic test

qrs <- function(meta,  mod) {
  mod.sig <- mod_sig.zR(meta, mod)
  mod.sig$ES <- (exp(2*mod.sig$ES)-1)/(exp(2*mod.sig$ES) + 1)
  mod.sig$Var <- (exp(2*mod.sig$Var)-1)/(exp(2*mod.sig$Var) + 1)
  mod.sig$Std.Err <-  sqrt(mod.sig$Var)
  k <-  mod.sig[mod.sig$Mod=="Overall", "K"] #sum(!is.na(mod)) 
  df <-  k-1
  mod.sig$z.value <- mod.sig$ES/mod.sig$Std.Err
  mod.sig$p.value <- round(2*(1-pt(abs(mod.sig$z.value),  df)), 6)
  mod.sig$L.95ci <- (exp(2*mod.sig$L.95ci)-1)/(exp(2*mod.sig$L.95ci) + 1)
  mod.sig$U.95ci <- (exp(2*mod.sig$U.95ci)-1)/(exp(2*mod.sig$U.95ci) + 1)
  mod.sig$df.Q <- mod.sig$K-1
  mod.sig$p_homog <-  ifelse(mod.sig$df.Q==0,  1,  round(pchisq(mod.sig$Q, mod.sig$df.Q, lower=FALSE), 5)) #work on this:if(mod.sig$df.Q[mod.sig$df.Q==0]) 1.0000 else pchisq(mod.sig$Q, mod.sig$df.Q, lower=FALSE)
  mod.sig$I2 <- (mod.sig$Q-(mod.sig$K-1))/mod.sig$Q      #I-squared is prop of total 
  mod.sig$I2 <- ifelse(mod.sig$I2<0, 0, mod.sig$I2)        #variation in the est of tmt
  #mod.sig$I2 <- ifelse(mod.sig$p_homog==1, mod.sig$I2==0,mod.sig$I2)
  mod.sig$I2 <- paste(round(mod.sig$I2*100, 4),  "%",  sep="")
  mod.sig2 <- mod.sig[, c(1, 2, 3, 6, 7, 9, 10, 4, 8, 5, 11, 12, 13)] #effects due to heterogeneity rather than chance (p263)
  rownames(mod.sig2) <- NULL
  colnames(mod.sig2) <- c("Mod", "k", "ES", "CI.lower", "CI.upper",  "Z", 
  "p",  "var.ES", "SE",  "Q", "df", "p.hetero",  "I2")
  return(mod.sig2)
}

#Q-statistic function (random effects)

CatModrQ <-  function(meta,  mod) {
  # Computes random effect Q-statistic (homogeneity test) for single predictor categorical 
  # moderator analysis. Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  # Returns:
  #   Random effects moderator Q-statistic, Q-within & between, df(Qw & Qb), and 
  #   homogeneity p-value within & between levels.
  mod.sig <- qrs(meta, mod)
  k <- mod.sig$k[mod.sig$Mod=="Overall"]                         #number of studies
  levels <- length(mod.sig$Mod)-1
  Qb.df <- levels-1 
  Qw.df <- k-levels
  Q <- mod.sig$Q[mod.sig$Mod=="Overall"]         #overall heterogeneity Q
  Qw <- sum(mod.sig$Q[!mod.sig$Mod=="Overall"])  #overall within-group heterogeneity statistic
  Qw_p.value <- 1-pchisq(Qw, Qw.df)      
  Qb <- Q-Qw                                     #overall between-group heterogeneity
  Qb_p.value <- 1-pchisq(Qb, Qb.df)
  mod.Qstat <- data.frame(Qb, Qb.df, Qb_p.value)
  names(mod.Qstat) <- c("Qb", "df.b", "p.b")
  return(mod.Qstat)
}

## new function 2.20.10
# Integrated function with fixed and random es, Q, etc
CatMod <- function(meta, mod) {
  fixed <- CatModf(meta,mod)
  random <-CatModr(meta,mod)
  random$Q <- fixed$Q
  random$df<- fixed$df
  random$p.hetero <- fixed$p.hetero
  random$I2 <- fixed$I2
  Qf <- CatModfQ(meta,mod) # fixed effect Q
  Qr <- CatModrQ(meta,mod)
  out<- list(Fixed=fixed, Q.fixed= Qf, Random=random, Q.random= Qr)
  return(out)
}


# Function for planned comparisons between 2 levels of moderator (random effects)

CatCompr <- function(meta,  mod, x1, x2,  method= "post.hoc1") {
  # Directly compares 2 levels of a categorical moderator using a random effects model. 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for moderator analysis.    
  #   x1: One level of categorical moderator
  #   x2: Other level (comparison group) of same categorical moderator
  #   method: "post.hoc" assumes the comparision was not planned prior to conducting
  #           the meta-analysis. The other option "planned" assumes you have planned 
  #           a priori to compare these levels of the categorical moderator. 
  #           Default is "post.hoc". 
  # Returns:
  #   Random effects moderator means per group, mean difference (d), variance of difference,
  #   p-value, and 95% confidence intervals. 
  modsig <- mod_sig.zR(meta, mod)
  com1 <- levels(modsig$Mod)[x1]  # first level of moderator
  com2 <- levels(modsig$Mod)[x2]  # second level of moderator
  x1.es <- modsig[modsig$Mod==com1, "ES"]
  x2.es <- modsig[modsig$Mod==com2, "ES"]
  x1.var <- modsig[modsig$Mod==com1, "Var"]
  x2.var <- modsig[modsig$Mod==com2, "Var"]
  g <- (-1)*x1.es + 1*x2.es  # pg 288 (Cooper et al.,  2009)
  var <- (-1)^2*x1.var + (1)^2*x2.var
  df <- 2-1
  chi.sqr <- g^2/var
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method == "post.hoc1") {  # post-hoc comparison (Tukey HSD method)
    fit<-aov(meta$z~meta$mod,weights=meta$wi.tau)
    fit<-TukeyHSD(fit)
  }
  if(method == "post.hoc2") {  # post-hoc comparison (Scheffe method)
    z <- g/sqrt(var)
    z2 <- z^2
    levels <- length(levels(as.factor(mod)))
    df.post <- levels-1 
    p.value <- 1-pchisq(z2, df.post)
    L.95ci <- g-1.96*sqrt(var)
    U.95ci <- g + 1.96*sqrt(var)
    fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
      names(fit) <- c("diff", "var.diff",  
      "p", "CI.lower", "CI.lower" )
  }
  if (method == "planned") {  # planned comparison (a priori)
   df <- 2-1
   chi.sqr <- g^2/var
   p.value <- 1-pchisq(chi.sqr, df)
   L.95ci <- g-1.96*sqrt(var)
   U.95ci <- g + 1.96*sqrt(var)
   fit <- data.frame(g, var, p.value, L.95ci, U.95ci)
   names(fit) <- c("diff", "var.diff",  
   "p", "CI.lower", "CI.lower" )
  }
  return(fit) 
}

#integrated catcomp function outputting both fixed and random

CatComp <- function(meta, mod, x1=NULL, x2=NULL,  method="post.hoc1") {
  fixed <- CatCompf(meta,mod, x1, x2, method)
  random <-CatCompr(meta,mod, x1, x2, method)
  out<- list(Fixed=fixed, Random=random)
  return(out)
}



# multifactor cat mod analysis [IN PROGRESS]:
#add relevant columns and convert back to r--it will works fine!
#MFCatMod <- function(meta, mod1, mod2) {
#  m <- wifun(meta)
#  m$mod1 <- mod1
#  m$mod2 <- mod2
#  fixed <- ddply(m, c("mod1", "mod2"), summarise, sum.wi = sum(wi),
#           sum.wiTi = sum(wiTi), sum.wiTi2 = sum(wiTi2))
#  fixed$ES <- fixed$sum.wiTi/fixed$sum.wi
#  random <- ddply(m, c("mod1", "mod2"), summarise, sum.wi.tau = sum(wi.tau),
#           sum.wiTi.tau = sum(wiTi.tau), sum.wiTi2.tau = sum(wiTi2.tau))
#  random$ES <- random$sum.wiTi.tau/random$sum.wi.tau
#  out <- list(Fixed = fixed, Random = random)
#  return(out)
#}

##==== META-REGRESSION FUNCTIONS (for continuous & categorical moderators)====##

# Meta-regression functions that correct the standard errors in OLS regressions 
# (Cooper,  2009; pp. 289-290). 
# These functions flexibly allows for single or multivariate predictors in the meta-regression 
# with continuous,  categorical,  or both moderator types simultaneously.

MAreg1 <- function(meta, mod, method="random") {  # Single predictor meta-regression
  # Computes single predictor fixed or random effects meta-regression (continuous or categorical). 
  # Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Moderator variable used for meta-regression.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  # Returns:
  #   Fixed or random effects beta coefficients,  adjusted standard errors, adjusted t-value, 
  #   95% confidence intervals, and adjusted p-value.
  m <- meta 
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- ddply(m,  .(id),  summarize,  r = aggrs(r), n=mean(n), mod = pick_one(mod))
  #meta <- agg_r2(m,m$mod)
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z <- 0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-(k - 1))/comp  # random effects variance
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method == "fixed") {
    reg0 <- lm(meta$z~1, weights=meta$wi)  # empty model
    reg <- lm(meta$z~mod, weights=meta$wi)  # model with mod  
    df <- anova(reg)["Residuals",  "Df"]
    ms.error <- anova(reg)["Residuals",  "Mean Sq"]
    t.crit <- qt(.975,  df)
    newSE.z <- round(summary(reg)$coef[, 2]/sqrt(ms.error), 6)  # 15.20 
    Bs.z <- round(summary(reg)$coef[, 1], 6)
    Bs <-  (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
    newSE <- (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
    t.adj <- round(Bs/newSE, 6)
    p.adj <- round(2*(1-pt(abs(t.adj),  df)), 6)
    lower.ci <- round(Bs-(t.crit*newSE), 6)  # 95% CI
    upper.ci <- round(Bs + (t.crit*newSE), 6)  # 95% CI
    modelfit <- MRfit(reg0,reg)   # internal function to assess fit
    #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
    #              cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
    #              symbols <- c("***",  "**",  "*",  ".",  " ")) 
  }
  # random effects #
  if(method == "random") {
    reg0 <- lm(meta$z~1, weights=meta$wi.tau) 
    reg <- lm(meta$z~mod, weights=meta$wi.tau)
    df  <- anova(reg)["Residuals",  "Df"]
    ms.error <- anova(reg)["Residuals",  "Mean Sq"]
    t.crit <- qt(.975,  df)
    newSE.z  <- round(summary(reg)$coef[, 2]/sqrt(ms.error), 6)  # 15.20 
    Bs.z  <- round(summary(reg)$coef[, 1], 6)  # Need an option to keep in z' 
    Bs <- (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
    newSE <- (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
    t.adj <- round(Bs/newSE, 6)
    p.adj <- round(2*(1-pt(abs(t.adj),  df)), 6)
    lower.ci <- round(Bs-(t.crit*newSE), 6)  # 95% CI
    upper.ci <- round(Bs + (t.crit*newSE), 6)  # 95% CI
    modelfit <- MRfit(reg0,reg)
    #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
    #              cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
    #              symbols <- c("***",  "**",  "*",  ".",  " ")) 
  }
  out <- data.frame(b=Bs, SE=newSE,  t=t.adj, CI.lower=lower.ci, CI.Upper=upper.ci, 
                      p=p.adj)
      
  return(list(out, modelfit))
}

##=== Multivariate Meta-Regression ===##

MAreg2  <-  function(reg) {  # Multivariate meta-regression
  # Computes multiple predictor fixed or random effects meta-regression (continuous and/or
  # categorical). Computations derived from chapter 15, Cooper et al. (2009). 
  # Args:
  #   reg: Weighted linear regression saved as an object (e.g., 
  #        reg <- lm(data$z ~ data$mod8 + data$mod1, weights= data$wi.tau). The outcome 
  #        variable is Fisher's z and the predictor moderators can be either continuous or
  #        categorical. Weight the regression by either the fixed or random effect weight
  #        (e.g., fixed=data$wi and random=data$wi.tau)
  # Returns:
  #   Fixed or random effects multivariate beta coefficients, adjusted standard errors, 
  #   adjusted t-value, 95% confidence intervals, and adjusted p-value.
  df  <-  anova(reg)["Residuals",  "Df"]
  ms.error  <-  anova(reg)["Residuals",  "Mean Sq"]
  t.crit  <-  qt(.975,  df)
  newSE.z  <-  summary(reg)$coef[, 2]/sqrt(ms.error)  
  Bs.z  <-  round(summary(reg)$coef[, 1], 4)
  Bs <-  (exp(2*Bs.z)-1)/(exp(2*Bs.z) + 1)
  newSE <-  (exp(2*newSE.z)-1)/(exp(2*newSE.z) + 1)
  t.adj  <-  round(Bs/newSE, 4)
  p.adj  <-  round(2*(1-pt(abs(t.adj),  df)), 4)
  lower.ci <- round(Bs-(t.crit*newSE), 4)     #95% CI
  upper.ci <- round(Bs + (t.crit*newSE), 4)     #95% CI
  #sig <- symnum(p.adj,  corr = FALSE,  na = FALSE,  
  #           cutpoints <- c(0,  0.001,  0.01,  0.05,  0.1,  1), 
  #           symbols <- c("***",  "**",  "*",  ".",  " ")) 
  out  <-  data.frame(b=Bs, SE=newSE,  t=t.adj, CI.lower=lower.ci, CI.Upper=upper.ci, 
                      p=p.adj)  
  return(out)
}


##============= GRAPHICS =============##

# requires ggplot2

##=== Meta-regression scatterplot with weighted regression line ===##

MAregGraph <- function(meta, mod,  method="random",  modname=NULL,  title=NULL, ylim=c(0, 1)) {
  # Outputs a scatterplot from a fixed or random effects meta-regression (continuous and/or
  # categorical). Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Moderator variable used for meta-regression.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  #   ylim: Limits of y-axis with the first argrument minimal value and second maximum value.
  #         Default is c(0,1).
  # Returns:
  #   Scatterplot with fixed or random effects regression line and where size of points are 
  #   based on study weights--more precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- meta$mod
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method=="fixed") {
    congraph <- ggplot(meta,  aes(mod, z, weight=wi), na.rm=TRUE) + 
    geom_point(aes(size=wi), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1), method = lm,  se = FALSE) + 
    xlab(modname) + ylab("Effect Size") +  
    ylim(ylim) +  
    opts(title=title, legend.position = "none")
  }
  if(method=="random") {
    congraph <- ggplot(meta,  aes(mod, z), na.rm=TRUE) + 
    geom_point(aes(size=wi.tau), alpha=.6, na.rm=TRUE) + 
    geom_smooth(aes(group=1, weight=wi.tau), method = lm, se = FALSE,  na.rm=TRUE) + 
    xlab(modname) + 
    ylab("Effect Size")  + 
    #ylim(min(meta$z), 1) +  
    opts(title=title, legend.position = "none")
  }
  return(congraph)
}

##=== Categorical Moderator Graph ===##

# Intermediate level function to add mean to boxplot

stat_sum_single1  <-  function(fun,  geom="point",  weight=wi, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",
                               geom=geom,  size = 7,  ...)      
}
stat_sum_single2  <-  function(fun,  geom="point",  weight=wi.tau, ...) {
                               stat_summary(fun.y=fun,  shape=" + ",  
                               geom=geom,  size = 7,  ...)      
}

CatModGraph <- function(meta, mod,  method="random",  modname=NULL,  title=NULL) {
  # Outputs a boxplot from a fixed or random effects moderator analysis.
  # Computations derived from chapter 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod: Categorical moderator variable used for analysis.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   modname: Name of moderator to appear on x-axis of plot. Default is NULL.
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Boxplot graph with median, interquartile range, max, min, and 
  #   outliers from a fixed or random effects categorical moderator analysis. Places
  #   jitter points for each study and the size of points are based on study weights--more 
  #   precise studies are larger. The ggplot2 package outputs the 
  #   rich graphics. 
  m <- meta
  m$mod <- mod
  compl<-!is.na(m$mod)
  meta<-m[compl,]
  #meta <- agg_r2(m,m$mod)
  #meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  #meta <- do.call(rbind, lapply(split(meta, meta$id), 
  #        function(.data) .data[sample(nrow(.data), 1),]))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  #meta <- meta
  #meta$mod <- mod
  #meta <- meta[!is.na(mod), ] 
  if(method=="fixed") {
    catmod <- ggplot(meta,  aes(factor(mod),  z, weight=wi), na.rm=TRUE) + 
                    geom_boxplot(outlier.size=2, na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) + 
                    stat_sum_single2(mean)

  }  
  if(method=="random") {
    catmod <- ggplot(meta,  aes(factor(mod),  z, weight=wi.tau), na.rm=TRUE) + 
                    geom_boxplot(outlier.size=2, na.rm=TRUE) + 
                    geom_jitter(aes(shape=factor(mod), size=wi.tau), alpha=.25) + 
                    #theme_bw() + 
                    xlab(modname) + 
                    ylab("Effect Size")  + 
                    opts(title=title, legend.position="none", na.rm=TRUE) + 
                    stat_sum_single2(mean)
  }
  return(catmod)
}

##=== Forrest Plot ===##

ForestPlot <- function(meta, method="random", title=NULL) {
  # Outputs a forest plot from a fixed or random effects omnibus analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size) for each study.
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Forest plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on 
  #   sample size) and 95% confidence intervals. The ggplot2 package outputs the rich graphics.  
  meta <- meta 
  meta$id <- factor(meta$id)  # , levels=rev(id))
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wi2 <- (1/meta$var.z)^2
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  sum.wi <- sum(meta$wi, na.rm=TRUE)        
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE) 
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)   
   if(method=="fixed") {  
     Tz.agg <- sum.wiTi/sum.wi              
     var.Tz.agg <- 1/sum.wi                  
     se.Tz.agg <- sqrt(var.Tz.agg)          
     T.agg <- (exp(2*Tz.agg)-1)/(exp(2*Tz.agg) + 1)  
     var.T.agg <- (exp(2*var.Tz.agg)-1)/(exp(2*var.Tz.agg) + 1) 
     Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                        
     k <- sum(!is.na(meta$r))                                   
     df <- k-1  
     omnibus <- data.frame(id="Omnibus", r=T.agg)
     meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create fixed ci for each study
     meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
     meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
     meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)
     forest <- ggplot(meta,  aes(y = id, x=r))+    #aes(y = factor(id, levels=rev(levels(id))),  x = r))  +  
                    geom_vline(xintercept=0) + 
                    geom_point(data=omnibus, colour="red", size=8, shape=23) + 
                    geom_point(aes(size=wi)) + 
                    opts(title=title,  legend.position="none") + 
                    geom_errorbarh(aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
                    geom_vline(colour="red", linetype=2,  xintercept=T.agg) + 
                    xlim(-1, 1) + 
                    xlab("Effect Size") + 
                    #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
                    ylab(NULL) 
  }
  if(method == "random") {
    comp <- sum.wi-sum.wi2/sum.wi
    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
    k <- sum(!is.na(meta$r))                                  
    df <- k-1      
    tau <- (Q-k + 1)/comp 				#random effects variance
    meta$var.tau <- meta$var.z + tau
    meta$wi.tau <- 1/meta$var.tau
    meta$wiTi.tau <- meta$wi.tau*meta$z
    meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE) 
    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
    Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
    var.Tz.agg.tau <-  1/sum.wi.tau                         #the following is inaccurate 14.23:  (Q - df)/comp
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    T.agg.tau <- (exp(2*Tz.agg.tau)-1)/(exp(2*Tz.agg.tau) + 1)
    var.T.agg.tau <- (exp(2*var.Tz.agg.tau)-1)/(exp(2*var.Tz.agg.tau) + 1)
    omnibus.tau  <- data.frame(id="Omnibus", r=T.agg.tau)
    meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
    meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
    meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
    meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)
    forest <- ggplot(meta,  aes(y = id,x = r))+    #aes(y = factor(id, levels=rev(levels(id))),  x = r))  +  
                  geom_vline(xintercept=0) + 
                  geom_point(data=omnibus.tau, colour="red",  size=8,  shape=23) + 
                  geom_point(aes(size=wi.tau)) + 
                  opts(title=title,  legend.position="none") + 
                  geom_errorbarh(aes(xmin = l.ci95,  xmax=u.ci95), size=.3, alpha=.6) + 
                  geom_vline(colour="red", linetype=2,  xintercept=T.agg.tau) + 
                  xlim(-1, 1) + 
                  xlab("Effect Size") + 
                  #scale_y_discrete(breaks = NA, labels=NA)+  # supress y-labels
                  ylab(NULL) 
  }
  return(forest)
}

##=== Funnel Plot ===## 

FunnelPlot <- function(meta, method="random",  title=NULL) {
  # Outputs a funnel plot from a fixed or random effects omnibus analysis to assess for
  # publication bias in the meta-analysis.
  # Computations derived from chapter 14, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size) for each study.
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Funnel plot with omnibus effect size (fixed or random), point for each study 
  #   where size of point is based on the study's precision (based primarily on sample 
  #   size) and standard error lines to assess for publication bias. 
  #   The ggplot2 package outputs the rich graphics.     
  meta <- MetaR(meta)
  sum.wi <- sum(meta$wi, na.rm=TRUE)       
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)   
  if(method=="fixed") {  
    meta$se.z <- sqrt(meta$var.z)           
    Tz.agg <- sum.wiTi/sum.wi                
    var.Tz.agg <- 1/sum.wi                  
    se.Tz.agg <- sqrt(var.Tz.agg)          
    Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                   
    k <- sum(!is.na(meta$r))                                  
    df <- k-1  
    l.ci95z <- meta$z-1.96*sqrt(meta$var.z)    
    u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
    omnibus.z <- Tz.agg
    funnel <- ggplot(meta,  aes(y = se.z,  x = z))  +  
                    geom_vline(colour="black", linetype=1, 
                    xintercept=omnibus.z) + 
                    geom_point(aes(size=wi)) + 
                    opts(title=title,  legend.position="none") + 
                    xlim(-1.7, 1.7) + 
                    ylim(.028, .5) + 
                    xlab("Fisher's z") + 
                    ylab("Standard Error") + 
                    stat_abline(intercept=omnibus.z/1.96, slope=(-1/1.96)) + 
                    stat_abline(intercept=(-omnibus.z/1.96), slope=1/1.96) + 
                    scale_y_continuous(trans="reverse")

  }
  if(method == "random") {
    meta$se.z.tau <- sqrt(meta$var.tau)
    sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)
    comp <- sum.wi-sum.wi2/sum.wi
    sum.wi.tau <- sum(meta$wi.tau, na.rm=TRUE)
    sum.wiTi.tau <- sum(meta$wiTi.tau, na.rm=TRUE)
    sum.wiTi2.tau <- sum(meta$wiTi2.tau, na.rm=TRUE)
    Tz.agg.tau <- sum.wiTi.tau/sum.wi.tau
    var.Tz.agg.tau <-  1/sum.wi.tau                        
    se.Tz.agg.tau <- sqrt(var.Tz.agg.tau)
    omnibus.z.tau  <- Tz.agg.tau
    l.ci95z <- meta$z-1.96*sqrt(meta$var.tau)    
    u.ci95z <- meta$z + 1.96*sqrt(meta$var.tau)
    l.ci95 <- (exp(2*l.ci95z)-1)/(exp(2*l.ci95z) + 1)
    u.ci95 <- (exp(2*u.ci95z)-1)/(exp(2*u.ci95z) + 1)
    meta$l.ci95 <- ifelse(l.ci95<=-1, -1, l.ci95)
    meta$u.ci95 <- ifelse(u.ci95>=1, 1, u.ci95)
    funnel <- ggplot(meta,  aes(y = se.z.tau,  x = z))  +  
                    geom_vline(colour="black", linetype=1, 
                    xintercept=omnibus.z.tau) + 
                    geom_point(aes(size=wi.tau)) + 
                    opts(title=title,  legend.position="none") + 
                    xlim(-1.7, 1.7) + 
                    ylim(.028, .5) + 
                    xlab("Fisher's z") + 
                    ylab("Standard Error") + 
                    stat_abline(intercept=omnibus.z.tau/1.96, slope=(-1/1.96)) + 
                    stat_abline(intercept=(-omnibus.z.tau/1.96), slope=1/1.96) + 
                    scale_y_continuous(trans="reverse")
  }
  return(funnel)
}

##== Multivariate Moderator Graphs ==##

MultiModGraph <- function(meta, conmod,  catmod, method="random",  
                          conmod.name=NULL,  title=NULL) {
  # Outputs a scatterplot and boxplot faceted by the categorical moderator from a 
  # fixed or random effects moderator analysis. Computations derived from chapter 
  # 14 and 15, Cooper et al. (2009). 
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   conmod: Continuous moderator variable used for analysis.
  #   catmod: Categorical moderator variable used for analysis.    
  #   method: Model used, either "random" or "fixed" effects. Default is "random".
  #   conmod.name: Name of continuous moderator variable to appear on x-axis of plot. 
  #   Default is NULL.
  #   title: Plot title. Default is NULL.
  # Returns:
  #   Multivariate moderator scatterplot graph faceted by categorical moderator levels. Also
  #   places a weighted regression line based on either a fixed or random effects analysis. 
  #   The ggplot2 packages outputs the rich graphics.  
  m <- meta
  m$conmod <- conmod
  m$catmod <- catmod
  compl <- !is.na(m$conmod)& !is.na(m$catmod)
  meta <- m[compl, ]
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  if(method=="fixed") {
    multimod <- ggplot(meta, aes(conmod, z, weight=wi), na.rm=TRUE) + 
                       opts(title=title, legend.position="none", na.rm=TRUE) + 
                       facet_wrap(~catmod)  + 
                       geom_point( aes(size=wi, shape=catmod)) + 
                       geom_smooth(aes(group=1, weight=wi),
                                   method= lm, se=FALSE, na.rm=TRUE) +
                       ylab("Effect Size") + 
                       xlab(conmod.name)
  }  
  if(method=="random") {
    multimod <- ggplot(meta, aes(conmod, z, weight=wi.tau), na.rm=TRUE) + 
                              opts(title=title, legend.position="none", na.rm=TRUE) + 
                              facet_wrap(~catmod)  + 
                              geom_point(aes(size=wi.tau, shape=catmod)) + 
                              geom_smooth(aes(group=1, weight=wi.tau), 
                                          method = lm, se = FALSE,  na.rm=TRUE) + 
                              ylab("Effect Size") + 
                              xlab(conmod.name)
  }
  return(multimod)
}

##===================== PUBLICATION BIAS ===================##
# Three approaches to assess for publication bias: (1) Fail
# Safe N, (2) Trim & Fill, and (3) Selection Modeling.

# Fail Safe N provides an estimate of the number of missing studies that 
# would need to exist to overturn the current conclusions.

PubBias <- function(meta, method = "fail.safe") {
  if(method == "fail.safe") {
    k <- length(meta$z.score)
    sum.z <- sum(meta$z.score)
    Z <- sum.z/sqrt(k)
    k0 <- round(-k + sum.z^2/(1.96)^2, 0)
    n.miss.per <- round(k0/k, 0)
  }
  out <- list("Fail.safe"=cbind(k,Z,k0, k.per))
  return(out)
}

##===================== PUBLICATION BIAS ===================##
# Three approaches to assess for publication bias: (1) Fail
# Safe N, (2) Trim & Fill, and (3) Selection Modeling.

# Fail Safe N provides an estimate of the number of missing studies that 
# would need to exist to overturn the current conclusions.

PubBias <- function(meta) {
  k <- length(meta$z.score)
  sum.z <- sum(meta$z.score)
  Z <- sum.z/sqrt(k)
  k0 <- round(-k + sum.z^2/(1.96)^2, 0)
  k.per <- round(k0/k, 0)
  out <- list("Fail.safe"=cbind(k,Z,k0, k.per))
  return(out)
}

##================== INTERRATER RELIABILITY ================##

# Kappa coefficients for inter-rater reliability (categorical variables)
# Imputs required are rater1 (first rater on Xi categorical variable)
# and rater2 (second rater on same Xi categorical variable)

Kappa <- function(rater1, rater2)  {
  # Computes Kappa coefficients for inter-rater reliability (categorical variables).
  # Args:
  #   rater1: First rater of categorical variable to be analyzed.
  #   rater2: Second rater on same categorical variable to be analyzed.
  # Returns:
  #   Kappa coefficients for inter-rater reliability (categorical variables).
  freq <- table(rater1, rater2)  # frequency table
  marg <- margin.table(freq) # total observations 
  marg2 <- margin.table(freq, 1)  # A frequencies (summed over rater2) 
  marg1 <- margin.table(freq, 2)  # B frequencies (summed over rater1)
  cellper <- prop.table(freq)  # cell percentages
  rowper <- margin.table(freq, 2)/margin.table(freq)  # row percentages 
  colper <- margin.table(freq, 1)/margin.table(freq)  # column percentages
  expected  <-  as.array(rowper) %*% t(as.array(colper)) 
  p.e <- sum(diag(expected))
  p.a_p.e <- sum(diag(cellper))- p.e
  p.e1 <- 1-p.e
  kappa <- p.a_p.e/p.e1
  return(kappa)
}

#Interclass correlations (ICC) for computing reliabilities for continuous variables
# Requires Psych package?

##====== Additional Functions ========##

# Function to reduce data set with complete data for 1 predictor 

ComplData <- function(meta, mod, type= "independent") {   
  # Outputs an aggregated data.frame that will remove any missing data from the data 
  # set. This is particularly useful to output non-missing data based on a specified
  # number of variables (generally in conjunction with the multivariate moderator
  # functions above)
  # Args:
  #   meta: data.frame with r (correlation coefficients) and n (sample size)for each study.
  #   mod1: Moderator variable wanting to be kept for further analysis.    
  # Returns:
  #   Reduced data.frame (with complete data) for the moderator entered into the 
  #   function while aggregating based on recommended procedured (Hunter & Schmidt,
  #   2004).  
  m <- meta
  m$mod <- mod
  compl <- !is.na(m$mod)
  m <- m[compl, ]
  #meta <- agg_r2(m,m$mod)
  if(type == "independent") {
    meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
    meta <- do.call(rbind, lapply(split(meta, meta$id), 
          function(.data) .data[sample(nrow(.data), 1),]))
  }
  if(type == "dependent") {
    meta <- ddply(m,  c("id", "mod"),  summarize,  r = aggrs(r), n=mean(n))
  }
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))   
  meta$var.z <- 1/(meta$n-3)
  meta$wi <-  1/meta$var.z
  meta$wiTi <- meta$wi*meta$z
  meta$wiTi2 <- meta$wi*(meta$z)^2
  meta$k <- length(meta$mod)
  mod <- meta$mod
  sum.wi <- sum(meta$wi, na.rm=TRUE)    
  sum.wi2 <- sum(meta$wi2, na.rm=TRUE)
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)    
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)
  comp <- sum.wi-sum.wi2/sum.wi
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi                         
  k <- sum(!is.na(meta$r))                                  
  df <- k-1      
  tau <- (Q-k + 1)/comp 				
  meta$var.tau <- meta$var.z + tau
  meta$wi.tau <- 1/meta$var.tau
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  return(meta)
}
         

# Correction for Attenuation

Rho_TU<- function(r,xx,yy) {
  r.corrected<-r/(sqrt(xx)*sqrt(yy))
  return(r.corrected)
   }


CorAtten <- function (meta, xx, yy) 
{
    m <- meta
    meta$xx <- xx
    meta$yy <- yy
    meta$r <- ifelse(is.na(meta$xx & meta$yy), meta$r, Rho_TU(meta$r, 
        meta$xx, meta$yy))
    meta$z <- 0.5 * log((1 + meta$r)/(1 - meta$r))
    meta$var.r <- ((1 - meta$r^2)^2)/(meta$n - 1)
    meta$var.r2 <- ifelse(is.na(m$xx & m$yy), meta$var.r, meta$var.r/(xx * 
        yy))
    meta$var.z <- 1/(meta$n - 3)
    meta$wi <- 1/meta$var.z
    meta$wiTi <- meta$wi * meta$z
    meta$wiTi2 <- meta$wi * (meta$z)^2
    meta$k <- length(meta$mod)
    mod <- meta$mod
    sum.wi <- sum(meta$wi, na.rm = TRUE)
    sum.wi2 <- sum(meta$wi2, na.rm = TRUE)
    sum.wiTi <- sum(meta$wiTi, na.rm = TRUE)
    sum.wiTi2 <- sum(meta$wiTi2, na.rm = TRUE)
    comp <- sum.wi - sum.wi2/sum.wi
    Q <- sum.wiTi2 - (sum.wiTi^2)/sum.wi
    k <- sum(!is.na(meta$r))
    df <- k - 1
    tau <- (Q - k + 1)/comp
    meta$var.tau <- meta$var.z + tau
    meta$wi.tau <- 1/meta$var.tau
    meta$wiTi.tau <- meta$wi.tau * meta$z
    meta$wiTi2.tau <- meta$wi.tau * (meta$z)^2
    return(meta)
}

# internal function for multifactor cat mod:
wifun <-  function(meta) {  
  meta$var.r <- var_r(meta$r, meta$n)
  meta$var.r <-  ifelse(is.na(meta$var.r), ((1-meta$r^2)^2)/(meta$n-1), meta$var.r)
  meta$z  <-  0.5*log((1 + meta$r)/(1-meta$r))  #computing r to z' for each study
  meta$var.z <- 1/(meta$n-3)  # computing var.z for each study
  meta$l.ci95z <- meta$z-1.96*sqrt(meta$var.z)     #create random ci for each study
  meta$u.ci95z <- meta$z + 1.96*sqrt(meta$var.z)
  meta$l.ci95 <- (exp(2*meta$l.ci95z)-1)/(exp(2*meta$l.ci95z) + 1)
  meta$u.ci95 <- (exp(2*meta$u.ci95z)-1)/(exp(2*meta$u.ci95z) + 1)  
  meta$z.score <- meta$z/sqrt(meta$var.z)
  meta$p.value <- 2*(1-pt(abs(meta$z.score),  meta$n-1))
  meta$wi <-  1/meta$var.z  # computing weight for each study
  meta$wiTi <- meta$wi*meta$z  # used to calculate omnibus
  meta$wiTi2 <- meta$wi*(meta$z)^2  # used to calculate omnibus
  # random effects #
  sum.wi <- sum(meta$wi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi <- sum(meta$wiTi, na.rm=TRUE)  # used to calculate random effects
  sum.wiTi2 <- sum(meta$wiTi2, na.rm=TRUE)  # used to calculate random effects
  Q <- sum.wiTi2-(sum.wiTi^2)/sum.wi  #used to calculate random effects            
  k <- sum(!is.na(meta$r))  # number of studies
  df <- k-1  # degree of freedom
  sum.wi2 <- sum(meta$wi^2, na.rm=TRUE)  # used to calculate random effects 
  comp <- sum.wi-sum.wi2/sum.wi  # (pg. 271) used to calculate random effects	
  meta$tau <- (Q-(k - 1))/comp  # Level 2 variance
  meta$var.tau <- meta$tau + meta$var.z  # Random effects variance (within study var + between var) 
  meta$wi.tau <- 1/meta$var.tau  # Random effects weights
  meta$wiTi.tau <- meta$wi.tau*meta$z
  meta$wiTi2.tau <- meta$wi.tau*(meta$z)^2
  meta$l.ci95z <- NULL 
  meta$u.ci95z <- NULL
 return(meta)
}

