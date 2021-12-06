\name{ProFound-package}
\alias{ProFound-package}
\alias{ProFound}
\docType{package}
\title{
\packageTitle{ProFound}
}
\description{
\packageDescription{ProFound}
}
\details{
\tabular{ll}{
Package: \tab ProFound\cr
Type: \tab Package\cr
Version: \tab 1.18.1\cr
Date: \tab 2021-12-06\cr
License: \tab LGPL-3\cr
Depends: \tab R (>= 3.1), FITSio, magicaxis (>= 2.0.8), Rcpp (>= 1.0.2)\cr
Imports: \tab data.table, celestial (>= 1.4.1), foreach\cr
Suggests: \tab ProFit, knitr, rmarkdown, EBImage, imager, LaplacesDemon, Rfast, fastmatch, snow, doSNOW, bigmemory, mvtnorm, Rwcs, Highlander (>= 0.1.7)\cr
}
}

\author{
\packageAuthor{ProFound}

Maintainer: \packageMaintainer{ProFound}
}
\references{
Robotham A.S.G., et al., 2018, MNRAS, 476, 3137
}
\keyword{profile}

\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
}
}
