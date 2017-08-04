\name{profoundMakeSigma}
\alias{profoundMakeSigma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make a Sigma Map
}
\description{
A utility function to construct a \code{ProFit} legal sigma map that can be input to \code{profitSetupData}.
}
\usage{
profoundMakeSigma(image, objects=0, sky=0, skyRMS=1, skycut=0, gain=1, readRMS=0, darkRMS=0,
image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU',
output_units='ADU', plot=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. Pixels set to 0 are interpreted as sky, and set to zero for calculating object shot-noise. If provided, this matrix *must* be the same dimensions as \option{image}. 
}
  \item{sky}{
Numeric; the absolute sky level. Consider using the sky output from \code{\link{profoundSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
Numeric; the RMS of the sky. Consider using the skyRMS output from \code{\link{profoundSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skycut}{
How many multiples of \option{skyRMS} above the \option{sky} to start calculating shot-noise based on the \option{gain} scaling of the \option{image}. If you are missing an object mask You almost certainly do not want this to be below 0 (else you will reduce the level of the sigma map just due to fluctuations in the sky), and in practice this should probably be set in the range 1-3.
}
  \item{gain}{
Numeric; the gain (in photo-electrons per ADU). For a very rough estimate consider using the gain output from \code{\link{profoundGainEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{readRMS}{
Numeric; the RMS of the read-noise. If you have estimated the sky RMS from the image directly this should not be necessary since it naturally captures this component. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{darkRMS}{
Numeric; the RMS of the dark-current-noise. If you have estimated the sky RMS from the image directly this should not be necessary since it naturally captures this component. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{image_units}{
Character; the units of the \option{image}. Must either be 'ADU' for generic astronomical data units, or 'elec' for photo-electrons.
}
  \item{sky_units}{
Character; the units of \option{sky} and \option{skyRMS}. Must either be 'ADU' for generic astronomical data units (the same type and scaling as per \option{image}), or 'elec' for photo-electrons.
}
  \item{read_units}{
Character; the units of \option{read}. Must either be 'ADU' for generic astronomical data units (the same type and scaling as per \option{image}), or 'elec' for photo-electrons.
}
  \item{dark_units}{
Character; the units of \option{dark}. Must either be 'ADU' for generic astronomical data units (the same type and scaling as per \option{image}), or 'elec' for photo-electrons.
}
  \item{output_units}{
Character; the units of the output sigma map. Must either be 'ADU' for generic astronomical data units (the same type and scaling as per \option{image}), or 'elec' for photo-electrons.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This is a simple utility function, but useful for beginners if they are unsure of how the error terms should be propagated (in short: in quadrature).
}
\value{
Numeric matrix; a sigma map the same size as \option{image}. This should be appropriate for feeding into \code{profitSetupData}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSkyEst}}, \code{\link{profoundGainEst}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits',
package="ProFound"))$imDat
sigma = readFITS(system.file("extdata", 'KiDS/G266035sigma.fits',
package="ProFound"))$imDat
segim = profoundSegImWatershed(image)
sigma_est=profoundMakeSigma(image, sky=segim$sky, skyRMS=segim$skyRMS, gain=1.8e13)

#This gets us close to the demo sigma image:
magplot(sigma, sigma_est, pch='.')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sigma }% use one of  RShowDoc("KEYWORDS")
