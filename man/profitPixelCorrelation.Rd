\name{profoundPixelCorrelation}
\alias{profoundPixelCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pixel to pixel correlation statistics
}
\description{
Returns the x and y dimension pixel-to-pixel correlation (often called covariance) at various scales, optionally returning a diagnostic plot.
}
\usage{
profoundPixelCorrelation(image, objects, mask, sky = 0, skyRMS = 1,
lag = c(1:9, 1:9 * 10, 1:9 * 100, 1:9 * 1000, 1:9 * 10000), fft = TRUE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as image.
}
  \item{mask}{
Boolean matrix; optional, parts of the \option{image} to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.  
}
  \item{sky}{
Numeric; the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{skyRMS}{
Numeric; the RMS of the sky. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{lag}{
Interger verctor; the pixel lags to measure pixel-to-pixel correlation over the x and y dimensions.
}
  \item{fft}{
Logical; if TRUE the 2D FFT is computed and the modulus image matrix is returned to \option{fft} and the (\option{image}-\option{sky})/\option{skyRMS} is return to \option{image_sky}, if FALSE the \option{fft} and \option{image_sky} objects are returned as NULL. \option{object} and \option{mask} pixels are used to identify pixels to replace as described below.
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
}
\details{
All statistics are computed on (\option{image}-\option{sky})/\option{skyRMS}. If \option{fft}=TRUE this matrix is return to \option{image_sky}.

The function is useful to assessing a number of image attributes. For one things it tells you whether all spatial variance has been detected and removed at small scales as objects (e.g. using \code{\link{profoundProFound}}), or at larger scales as sky fluctuations. Assuming the object detection and sky removal has worked well, the remaining pixel-to-pixel correlation likely represents instrument level covariance. In practice nearly all processes produce positive pixel correlation, but it is not impossible that negative correlation can be introduced during the reduction process, particularly when over-subtracting the sky around bright stars.

For calculating the raw pixel-to-pixel correlation (as returned by \option{cortab}) \option{mask} and \option{object} pixels are ignored, so correlation is only considered where both pixels are flagged as un-masked sky pixels. The 2D image FFT output (option{fft}) replaces masked or object pixels with Normally distributed noise after the input \option{image} has had the \option{sky} subtracted and divided by the \option{skyRMS}. Note that this means the FFT generated is partly stochastic (it will differ a bit each time it is run), but in practice it will be quite persistant for small scales (the centre) and stochastic at larger scales (around the edge of the FFT image).

The relative standard-deviations returned in \option{cortab} are calculated by taking the standard-deviation of the lagged pixel differences of (\option{image}-\option{sky})/\option{skyRMS} and dividing through by sqrt(2). This means for well behaved data they should be 1, and the dashed lines on the diagnostic plot should fall on 1.
}
\value{
A list containing three objects:

\item{cortab}{A five column data.frame, containing the pixel lag (lag), the correlation in the x-dimension (corx), the correlation in the y-dimension (cory), the pixel lag implied relative standard-deviation in the x-dimension (relsdx) and the pixel lag implied relative standard-deviation in the y-dimension (relsdy).}
\item{fft}{Numeric matrix; if \option{fft}=TRUE this object contains the modulus of the 2D FFT of the \option{image} with the same dimensions, if \option{fft}=FALSE it is NULL. The output is centred such that small scale structures are in the centre of the image matrix.}
\item{image_sky}{Numeric matrix; if \option{fft}=TRUE this object contains the (\option{image}-\option{sky})/\option{skyRMS}, if \option{fft}=FALSE it is NULL.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundProFound}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

corout_raw=profoundPixelCorrelation(image$imDat, plot=TRUE)
magimage(corout_raw$fft)
points(178, 178, cex=10, col='red')

#There is clearly some residual structure masking out the brighter parts of objects:

corout_objects=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects, plot=TRUE)
magimage(corout_objects$fft)
points(178, 178, cex=10, col='red')

#Using the more aggressive objects_redo removed nearly all of this:

corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo, plot=TRUE)
magimage(corout_objects_redo$fft)
points(178, 178, cex=10, col='red')

#We can use the pixel correlation function, in particular the FFT output, to assess how
#much further we can afford to push the source extraction in our image.

profound=profoundProFound(image, skycut=2.0, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=1.0, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=0.8, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=0.6, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft, verbose=TRUE, plot=TRUE)

#By doing ProFoundsource detection on the FFT itself it tells us if there are significant
#sources of a certain common scale (usually small) still in the image to extract.
#The levels above suggest we cannot push much further than a skycut=1.0. Clearly using
#skycut=0.6 introduces a lot of fake sources.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ FFT }% __ONLY ONE__ keyword per line
