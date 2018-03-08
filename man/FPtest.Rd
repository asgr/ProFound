\name{FPtest}
\alias{FPtest}
\docType{data}
\title{
False Positive Reference Data
}
\description{
This data consists of 1,000 runs of a random 1000 x 1000 noise matrix through \code{\link{profoundProFound}}. The catalogue is a concatenation of all the segstats outputs for all of these run.
}
\usage{data("FPtest")}
\format{
  A data frame with 7012 observations on the following 56 variables. See \code{\link{profoundProFound}} for a detailed discussion on each of these parameters.
  \describe{
    \item{\code{segID}}{a numeric vector}
    \item{\code{uniqueID}}{a numeric vector}
    \item{\code{xcen}}{a numeric vector}
    \item{\code{ycen}}{a numeric vector}
    \item{\code{xmax}}{a numeric vector}
    \item{\code{ymax}}{a numeric vector}
    \item{\code{RAcen}}{a logical vector}
    \item{\code{Deccen}}{a logical vector}
    \item{\code{RAmax}}{a logical vector}
    \item{\code{Decmax}}{a logical vector}
    \item{\code{sep}}{a numeric vector}
    \item{\code{flux}}{a numeric vector}
    \item{\code{mag}}{a numeric vector}
    \item{\code{cenfrac}}{a numeric vector}
    \item{\code{N50}}{a numeric vector}
    \item{\code{N90}}{a numeric vector}
    \item{\code{N100}}{a numeric vector}
    \item{\code{R50}}{a numeric vector}
    \item{\code{R90}}{a numeric vector}
    \item{\code{R100}}{a numeric vector}
    \item{\code{SB_N50}}{a numeric vector}
    \item{\code{SB_N90}}{a numeric vector}
    \item{\code{SB_N100}}{a numeric vector}
    \item{\code{xsd}}{a numeric vector}
    \item{\code{ysd}}{a numeric vector}
    \item{\code{covxy}}{a numeric vector}
    \item{\code{corxy}}{a numeric vector}
    \item{\code{con}}{a numeric vector}
    \item{\code{asymm}}{a logical vector}
    \item{\code{flux_reflect}}{a logical vector}
    \item{\code{mag_reflect}}{a logical vector}
    \item{\code{semimaj}}{a numeric vector}
    \item{\code{semimin}}{a numeric vector}
    \item{\code{axrat}}{a numeric vector}
    \item{\code{ang}}{a numeric vector}
    \item{\code{signif}}{a numeric vector}
    \item{\code{FPlim}}{a numeric vector}
    \item{\code{flux_err}}{a numeric vector}
    \item{\code{mag_err}}{a numeric vector}
    \item{\code{flux_err_sky}}{a numeric vector}
    \item{\code{flux_err_skyRMS}}{a numeric vector}
    \item{\code{flux_err_shot}}{a numeric vector}
    \item{\code{sky_mean}}{a numeric vector}
    \item{\code{sky_sum}}{a numeric vector}
    \item{\code{skyRMS_mean}}{a numeric vector}
    \item{\code{Nedge}}{a logical vector}
    \item{\code{Nsky}}{a logical vector}
    \item{\code{Nobject}}{a logical vector}
    \item{\code{Nborder}}{a logical vector}
    \item{\code{Nmask}}{a logical vector}
    \item{\code{edge_frac}}{a logical vector}
    \item{\code{edge_excess}}{a logical vector}
    \item{\code{flag_border}}{a logical vector}
    \item{\code{iter}}{a numeric vector}
    \item{\code{origfrac}}{a numeric vector}
    \item{\code{flag_keep}}{a logical vector}
  }
}
\details{
Specifically we ran with defaults the following command 1,000 times in a loop:

profoundProFound(matrix(rnorm(1e6),1e3))

The output is then a reference of the false positive rate, since we have not injected any sources into the images. The fact we find 7,012 false detections mean we expect 7 false positives per 1e6 pixels (the size in pixels of the input matrix). To compare against any target data we need to adjust the magnitudes by the sky RMS magnitude level, i.e. add on profoundFlux2Mag(skyRMS, 0) (if the zero point is 0 for our target data). See Examples for a comparison to our included VIKING data.
}
\source{
FPtest={}

for(i in 1:1000){FPtest=rbind(FPtest,profoundProFound(matrix(rnorm(1e6),1e3))$segstats)}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, magzero=30, rotstats=TRUE)
skyRMS=median(profound$skyRMS)
magoff=profoundFlux2Mag(skyRMS, 30)
totpix=prod(profound$dim)

#We can easily compute the expected number of false positives on an image this size:
data("FPtest")
dim(FPtest)[1]*totpix/1e6/1e3

#And plot the detections and expected false positive distributions:
maghist(profound$segstats$mag, seq(-11,-1,by=0.2)+magoff)
maghist(FPtest$mag+magoff, seq(-6,-1,by=0.2)+magoff, scale=totpix/1e6/1e3, add=TRUE,
border='red')
}
}
\keyword{datasets}
