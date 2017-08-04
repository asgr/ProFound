\name{ProFound-package}
\alias{ProFound-package}
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
Version: \tab 1.0.0\cr
Date: \tab 2017-07-24\cr
License: \tab GPL-3+\cr
Depends: \tab R (>= 3.0), magicaxis (>= 2.0.3), celestial (>= 1.4.1), FITSio, data.table\cr
Suggests: \tab ProFit, knitr, rmarkdown, foreach, doParallel, EBImage, akima, imager\cr
}
}

\author{
\packageAuthor{ProFound}

Maintainer: \packageMaintainer{ProFound}
}
\references{
Robotham A.S.G., et al., 2017, MNRAS, 466, 1513
}
\keyword{profile}

\examples{

# Load ProFound example data

image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))

profound=profitProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

}
