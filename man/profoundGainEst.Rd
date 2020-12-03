\name{profoundGainEst}
\alias{profoundGainEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Gain Estimator
}
\description{
High level function to estimate a rough value for the image gain in cases where you have no idea what the true image gain is. In practice this tends to be accurate to an order of magnitude and provides a reasonable lower limit for the true gain, which is good enough to make a rough first attempt at a sigma map.
}
\usage{
profoundGainEst(image = NULL, mask = 0, objects = 0, sky = 0, skyRMS = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. Pixels set to 0 are interpreted as sky, and set to zero for calculating object shot-noise. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
Numeric scalar; user provided estimate of the absolute sky level. If this is not provided then it will be computed internally using \code{\link{profoundSkyEst}}.
}
  \item{skyRMS}{
Numeric scalar; user provided estimate of the RMS of the sky. If this is not provided then it will be computed internally using \code{\link{profoundSkyEst}}.  
}
}
\details{
This function makes use of the fact that a true Poisson distribution cannot generate samples below 0 and the distribution shape properties of the sky pixels. In practice this means the gain estimated is low as it can be. Once the ProFit fit has been made the gain estimated can be improved based on the residuals (assuming the model does a good job of subtracting the data).
}
\value{
Numeric scalar; the estimated gain of the \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundMakeSegim}}, \code{\link{profoundMakeSegimExpand}}, \code{\link{profoundSkyEst}}, \code{\link{profoundMakeSigma}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))
profound=profoundProFound(image)
profoundGainEst(image$imDat, objects=profound$objects_redo, sky=profound$sky)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ gain }% use one of  RShowDoc("KEYWORDS")
