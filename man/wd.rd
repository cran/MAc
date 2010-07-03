\name{wd}
\alias{wd}
\title{Output to Word in formatted tables 
}
\description{Function for exporting MA output to nicely formatted Word tables.
}
\usage{
wd(object, get = FALSE, new = FALSE, ...)
}
\arguments{
  \item{object}{Will take either an \code{omni} (Omnibus), \code{mareg} (meta-regression), or \code{macat} (single predictor categorical moderator analysis) object and export to Word in a formatted table. 
}
  \item{get}{Start up the Word program? TRUE if an instance of Word is not currently open.
}
  \item{new}{Output data into a new Word document? TRUE or FALSE. 
}
 \item{...}{ Additional arguments to be passed to R2wd functions.
  }   
}
\details{This function depends of the \code{R2wd} package. See Christian Ritter (2010). R2wd: Write MS-Word documents from R. R
  package version 1.3 for the details of the \code{R2wd} package.
}
\references{ Christian Ritter (2010). R2wd: Write MS-Word documents from R. R
  package version 1.3.
} 
\author{AC Del Re & William T. Hoyt

Maintainer: AC Del Re \email{acdelre@gmail.com}
}
\seealso{
\code{\link{omni}},
\code{\link{mareg}},
\code{\link{macat}}
}
\examples{
# Sample data
id<-c(1:20)
n.1<-c(10,20,13,22,28,12,12,36,19,12,36,75,33,121,37,14,40,16,14,20)
n.2 <- c(11,22,10,20,25,12,12,36,19,11,34,75,33,120,37,14,40,16,10,21)
g <- c(.68,.56,.23,.64,.49,-.04,1.49,1.33,.58,1.18,-.11,1.27,.26,.40,.49,
.51,.40,.34,.42,1.16)
var.g <- c(.08,.06,.03,.04,.09,.04,.009,.033,.0058,.018,.011,.027,.026,.0040,
.049,.0051,.040,.034,.0042,.016)
mod<-factor(c(rep(c(1,1,2,3),5)))
mods2<-c(rep(1:5,4))
df<-data.frame(id, n.1,n.2, g, var.g,mod, mods2)


# Examples

# install R2wd
# install.packages('R2wd', dependencies = TRUE)

# mareg fuction
temp <- mareg(g~ mod + mods2, var = var.g, method = "REML", data = df)

# Export data to Word in formatted table
\dontrun{wd(temp, get = TRUE, new = TRUE)}
}
\keyword{word}

