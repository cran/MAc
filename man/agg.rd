\name{agg}
\alias{agg}
\title{Meta-Analysis Aggregation
}
\description{This fuction will simultaneously aggregate all within-study effect sizes while taking into account the correlations among the within-study outcome measures (Hunter and Schmidt, 2004). The default imputed correlation between within-study effect sizes is set at .50 (Wampold et al., 1997) and will compute an aggregated effect size for each study. This default of .50 is adjustable.  
}
\usage{
agg(id, r, n, cor = .50, mod=NULL, data)
}
\arguments{
 \item{id}{Study id with multiple rows of same id.
}
\item{r}{Correlation coefficient.
}
\item{n}{Sample size.
}
\item{cor}{Estimated correlation among within-study outcome variables. Default is .50 based on Wampold et al. (1997) recommended procedure.
}
\item{mod}{Default is NULL. To aggregate by id and one moderator. If there are multiple levels of a categorical moderator within study and one can in derive seperate effect size estimates for each level within and between studies. However, there will be dependency issues and one way to resolve is shown below in the examples. 
}
\item{data}{\code{data.frame} with above values.
}
}
\value{Outputs a \code{data.frame} with aggregated effect sizes where each study is reduced to one row per study (unless aggregated by a moderator) by a weighted average formula. This formula is based on Hunter and Schmidt's (2004) approach to aggregation of dependent effect sizes (see chapter 10, pp. 435-8). 
}
\references{Hunter, J. E., Schmidt, F. L. (2004). Methods of meta-analysis (2nd edition). Thousand Oaks, CA: Sage.
}

\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\examples{
id<-rep(1:5, 4)
n<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
r<-c(.68,.56,.23,.64,.49,-.04,.49,.33,.58,.18,-.11,.27,.26,.40,.49,
.51,.40,.34,.42,.16)
mod1<-factor(rep(1:2, 10))
dat<-data.frame(id,n,r,mod1)

# Examples

# aggregate to 1 id per study (independent sample)
agg(id = id, r = r, n = n, data=dat) 

# aggregate by id & a moderator (non-independent sample)
temp <- agg(id = id, r = r, n = n, mod = mod1, data=dat) 

temp

# This function below will randomly select one within
# study level of the moderator (if there are more than one) and output an
# independent sample. Replace temp with the name of your data.
do.call(rbind, lapply(split(temp, temp$id), 
          function(.data) .data[sample(nrow(.data), 1),]))


}
\keyword{ aggregation }

