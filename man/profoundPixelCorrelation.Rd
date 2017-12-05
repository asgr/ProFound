\name{profoundPixelCorrelation}
\alias{profoundPixelCorrelation}
\alias{profoundSkySplitFFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pixel to pixel correlation statistics
}
\description{
Returns the x and y dimension pixel-to-pixel correlation (often called covariance) at various scales, optionally returning a diagnostic plot.
}
\usage{
profoundPixelCorrelation(image, objects, mask, sky = 0, skyRMS = 1,
lag = c(1:9, 1:9 * 10, 1:9 * 100, 1:9 * 1000, 1:9 * 10000), fft = TRUE, plot = FALSE,
ylim=c(-1,1), log='x', grid=TRUE, ...)
profoundSkySplitFFT(image, objects, mask, sky = 0, skyRMS = 1, skyscale=100)
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
Logical; should a x/y correlation diagnostic plot be generated?
}
  \item{ylim}{
Numeric vector; range of data to display (see \code{\link{magplot}} for details). Only relevant if \option{plot}=TRUE.
}
  \item{log}{
Character scalar; log axis arguments to be passed to plot. E.g. use 'x', 'y', 'xy' or 'yx' as appropriate (see \code{\link{magplot}} for details). Only relevant if \option{plot}=TRUE.
}
  \item{grid}{
Logical; indicates whether a background grid should be drawn onto the plotting area (see \code{\link{magplot}} for details). Only relevant if \option{plot}=TRUE.
}
  \item{skyscale}{
Numeric scalar; required, the pixel scale that the FFT should split the provided \option{image_sky} at. This should be chosen so as to separate out true sky modes and possible sources still in the sky. Too small and real sources will be put into the \option{sky_lo} image returned, so larger is usually safer.
}
  \item{\dots}{
Further arguments to parse to \code{\link{magplot}}. Only relevant if \option{plot}=TRUE.
}
}
\details{
\code{profoundPixelCorrelation}:

All statistics are computed on (\option{image}-\option{sky})/\option{skyRMS}. If \option{fft}=TRUE this matrix is return to \option{image_sky}.

The function is useful to assessing a number of image attributes. For one things it tells you whether all spatial variance has been detected and removed at small scales as objects (e.g. using \code{\link{profoundProFound}}), or at larger scales as sky fluctuations. Assuming the object detection and sky removal has worked well, the remaining pixel-to-pixel correlation likely represents instrument level covariance. In practice nearly all processes produce positive pixel correlation, but it is not impossible that negative correlation can be introduced during the reduction process, particularly when over-subtracting the sky around bright stars.

For calculating the raw pixel-to-pixel correlation (as returned by \option{cortab}) \option{mask} and \option{object} pixels are ignored, so correlation is only considered where both pixels are flagged as un-masked sky pixels. The 2D image FFT output (\option{fft}) replaces masked or object pixels with Normally distributed noise after the input \option{image} has had the \option{sky} subtracted and divided by the \option{skyRMS}. Note that this means the FFT generated is partly stochastic (it will differ a bit each time it is run), but in practice it will be quite persistant for small scales (the centre) and stochastic at larger scales (around the edge of the FFT image).

The relative standard-deviations returned in \option{cortab} are calculated by taking the standard-deviation of the lagged pixel differences of (\option{image}-\option{sky})/\option{skyRMS} and dividing through by sqrt(2). This means for well behaved data they should be 1, and the dashed lines on the diagnostic plot should fall on 1.

\code{profoundSkySplitFFT}:

The FFT split output separates the provided image into hi frequency (\option{sky_hi}) and low frequency (\option{sky_lo}) modes. The idea is that \option{sky_lo} might represent additional sky with complex structure (not captured by the bicubic/bilinear extimated sky) that still needs to be subtracted off the image, whilst \option{sky_hi} might contain some as yet un-subtracted sources.

In principle \code{profoundSkySplitFFT} can be run with any image, but the separation into the low and high frequency modes is not easily interpretable in the presence of many real objects since they will dominate the power at all scales (trust me on this).
}
\value{
\code{profoundPixelCorrelation}:

A list containing three objects:

\itemize{
  \item cortab: A data.frame containing: 
  \itemize{
    \item lag: The pixel lag
    \item corx: The correlation in the x-dimension
    \item cory: The correlation in the y-dimension
    \item corx_neg: The correlation of sub sky versus sky pixels in x
    \item cory_neg: The correlation of sub sky sversus ky pixels in y
    \item corx_pos: The correlation of excess sky versus sky pixels in x
    \item cory_pos: The correlation of excess sky versus sky pixels in y
    \item corx_diff: corx_pos - corx_neg
    \item cory_diff: cory_pos - cory_neg
    \item relsdx: The pixel lag implied relative standard-deviation in x
    \item relsdy: The pixel lag implied relative standard-deviation in y
  }
  \item fft: if \option{fft}=TRUE this object contains a list containing x, y, and z. If \option{fft}=FALSE it is NULL. x and y contain the frequency values of the 2D FFT in Hz (the DC or 0 Hz signal is in the absolute centre). z contains the power component of the 2D FFT image as a numeric matrix; the modulus of the 2D FFT of the \option{image} with the same dimensions. The output is centred such that large scale structures are in the centre of the image matrix.
  \item image_sky: Numeric matrix; if \option{fft}=TRUE this object contains the (\option{image}-\option{sky})/\option{skyRMS}, if \option{fft}=FALSE it is NULL.
  
\code{profoundSkySplitFFT}:

A list containing three numeric matrices:

\item{sky}{The new sky estimate, defined as the input \option{sky}+\option{sky_lo}.}
\item{sky_lo}{The low frequency modes extracted from the objects masked \option{image}-\option{sky}.}
\item{sky_hi}{The high frequency modes extracted from the objects masked \option{image}-\option{sky}.}
}
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
magimage(corout_raw$fft, xlab='x / Hz', ylab='y / Hz')
points(0, 0, cex=10, col='red')

#There is clearly some residual structure masking out the brighter parts of objects:

corout_objects=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects, plot=TRUE)
magimage(corout_objects$fft, xlab='x / Hz', ylab='y / Hz')
points(0, 0, cex=10, col='red')

#Using the more aggressive objects_redo removed nearly all of this:

corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo, plot=TRUE)
magimage(corout_objects_redo$fft, xlab='x / Hz', ylab='y / Hz')
points(0, 0, cex=10, col='red')

#We can use the pixel correlation function, in particular the FFT output, to assess how
#much further we can afford to push the source extraction in our image.

profound=profoundProFound(image, skycut=2.0, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft$z, skycut=2, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft$z, skycut=2, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=1.0, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft$z, skycut=2, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=0.8, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft$z, skycut=2, verbose=TRUE, plot=TRUE)

profound=profoundProFound(image, skycut=0.6, magzero=30, verbose=TRUE, plot=TRUE)
corout_objects_redo=profoundPixelCorrelation(image$imDat, sky=profound$sky,
skyRMS=profound$skyRMS, objects=profound$objects_redo)
magimage(corout_objects_redo$image_sky)
profoundProFound(corout_objects_redo$fft$z, skycut=2, verbose=TRUE, plot=TRUE)

#By doing ProFoundsource detection on the FFT itself it tells us if there are significant
#sources of a certain common scale (usually small) still in the image to extract.
#The levels above suggest we cannot push much further than a skycut=1.0. Clearly using
#skycut=0.6 introduces a lot of fake sources.

#We can improve the sky using profoundSkySplitFFT

profound=profoundProFound(image, type="bicubic")
newsky=profoundSkySplitFFT(image$imDat, objects=profound$objects_redo, sky=profound$sky,
skyRMS=profound$skyRMS)

#Old versus new sky:

magimage(profound$sky)
magimage(newsky$sky)

#Original image, old sky subtraction and new sky subtraction (pretty subtle!):

magimage(image$imDat)
magimage(image$imDat-profound$sky)
magimage(image$imDat-newsky$sky)

#Be warned, you need a reasonable estimate of the sky and objects before running this.
#If we run on the original image that even the high/low frequency modes look very odd:

magimage(profoundSkySplitFFT(image$imDat)$sky_lo)
magimage(profoundSkySplitFFT(image$imDat)$sky_hi)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }% use one of  RShowDoc("KEYWORDS")
\keyword{ FFT }% __ONLY ONE__ keyword per line
