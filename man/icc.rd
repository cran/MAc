\name{icc}
\alias{icc}
\title{Intraclass correlation coefficient (ICC) for oneway and twoway models
}
\description{Computes single score or average score ICCs as an index of interrater reliability of quantitative data. Additionally, F-test and confidence interval are computed.
}
\usage{icc(ratings, model = c("oneway", "twoway"), 
    type = c("consistency", "agreement"), 
    unit = c("single", "average"), r0 = 0, conf.level = 0.95)
}
\arguments{
\item{ratings}{	 n*m matrix or dataframe, n subjects m raters.
}
\item{model}{a character string specifying if a \code{oneway} model (default) with row effects random, or a \code{twoway} model with column and row effects random should be applied. You can specify just the initial letter.
}
\item{type}{	 a character string specifying if '"consistency"' (default) or '"agreement"' between raters should be estimated. If a '"oneway"' model is used, only '"consistency"' could be computed. You can specify just the initial letter.
}
\item{unit}{a character string specifying the unit of analysis: Must be one of \code{single} (default) or \code{average}. You can specify just the initial letter.
}
\item{r0}{	 specification of the null hypothesis r = r0. Note that a one sided test (H1: r > r0) is performed.
}
\item{conf.level}{confidence level of the interval.
}
}

\details{
This function was created by Matthias Gamer for the \code{irr} package. For more details, see: 

http://rss.acs.unt.edu/Rdoc/library/irr/html/icc.html

Details for the function:

Missing data are omitted in a listwise way. When considering which form of ICC is appropriate for an actual set of data, one has take several decisions (Shrout & Fleiss, 1979):

1. Should only the subjects be considered as random effects (\code{oneway} model) or are subjects and raters randomly chosen from a bigger pool of persons (\code{twoway} model).

2. If differences in judges' mean ratings are of interest, interrater \code{agreement} instead of \code{consistency} should be computed.

3. If the unit of analysis is a mean of several ratings, unit should be changed to \code{average}. In most cases, however, single values (\code{unit = single}) are regarded.
}
\value{A list with class \code{icclist} containing the following components:
\item{subjects}{	 the number of subjects examined.
}
\item{raters}{	 the number of raters.
}
\item{model}{	 a character string describing the selected model for the analysis.
}
\item{type}{	 a character string describing the selected type of interrater reliability.
}
\item{unit}{	 a character string describing the unit of analysis.
}
\item{icc.name}{	 a character string specifying the name of ICC according to McGraw & Wong (1996).
}
\item{value}{	 the intraclass correlation coefficient.
}
\item{r0}{	 the specified null hypothesis.
}
\item{Fvalue}{	 the value of the F-statistic.
}
\item{df1}{the numerator degrees of freedom.
}
\item{df2}{	 the denominator degrees of freedom.
}
\item{p.value}{	 the p-value for a two-sided test.
}
\item{conf.level}{	 the confidence level for the interval.
}
\item{lbound}{	 the lower bound of the confidence interval.
}
\item{ubound}{	 the upper bound of the confidence interval.
}
}
\author{Matthias Gamer
}
\references{

Bartko, J.J. (1966). The intraclass correlation coefficient as a measure of reliability. Psychological Reports, 19, 3-11.

McGraw, K.O., & Wong, S.P. (1996), Forming inferences about some intraclass correlation coefficients. Psychological Methods, 1, 30-46.

Shrout, P.E., & Fleiss, J.L. (1979), Intraclass correlation: uses in assessing rater reliability. Psychological Bulletin, 86, 420-428.
}
\examples{
# sample data

study <- c(1,1,2,2,3,3)
rater <- c(rep(1:2,3))
mod1 <- round(rnorm(6, 10, 1))
mod2 <- c(5,5, 9, 9, 8, 8)
mod3 <- c(10,10, 9, 9, 8, 8)
w <-data.frame(study, rater, mod1, mod2, mod3)
w


# if data is in this format:

# study rater mod1 mod2 mod3
#     1     1    9    9   10
#     1     2   11    8   10
#     2     1    9   10   11
#     2     2    9   10   11
#     3     1    9    9    8
#     3     2   12    9    8
#
# the data will need to be reshaped to be processed by the 
# icc function:

long <- reshape(w, varying=colnames(w)[3:5], v.names="Code", 
            idvar=c('study', 'rater'), timevar="mods", direction='long')
wide <- reshape(long, idvar=c('mods', 'study'), timevar='rater')

# icc function (created by Matthias Gamer for the 'irr' package)

icc(cbind(wide$Code.1, wide$Code.2), type= "consistency")
}
