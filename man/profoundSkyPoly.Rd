\name{profoundSkyPoly}
\alias{profoundSkyPoly}
\alias{profoundSkyPlane}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit N-order Polynomial Sky
}
\description{
Functions to fit polynomial sky planes to image data. Includes a useful summary plot. These might be preferable in complicated situations where a large amount of the image is masked by source/s, e.g. a large foreground cluster in the centre of the image. In these cases the restrictions of the fits here might interpolated over masked regions more nicely than the default \code{\link{profoundMakeSkyGrid}} behaviour. The downside is much more care is needed to mask out the image appropriately, hence the inclusion of a quancut argument and diagnostic plots.
}
\usage{
profoundSkyPoly(image, objects = NULL, degree = 1, quancut = NULL, mask = NULL, 
  mode_shift = FALSE, plot = FALSE, ...)
profoundSkyPlane(image, objects = NULL, quancut = NULL, mask = NULL, mode_shift = FALSE, 
  plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix (1,0); optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}. If provided then this is used to initially mask pixels for determining the correct \option{sky}.
}
  \item{degree}{
Integer scalar; the polynomial degree to model the sky with. \option{degree} = 0 is a pedestal sky, 1 is a linear plane, 2 is quadratic, 3 cubic, 4 quartic and 5 quintic etc. As a convenience setting \option{degree} = -1 sets the sky to be 0 throughout (perhaps useful for checking diagnostics).
}
  \item{quancut}{
Numeric scalar; optional, the level at which to mask out brighter pixels for analysis. This might be useful if \option{objects} has not already been computed.
}
  \item{mask}{
Boolean matrix or integer scalar (1,0); optional, parts of the image to mask out (i.e. ignore). If a matrix is provided, this matrix *must* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis.
}
  \item{mode_shift}{
Logical; should the mode be shifted to zero regardless of best polynomical regression found? As long as the non-\option{object}, non-\option{mask} and non-\option{quancut} excluded pixels contain the true sky mode, this is usually a good idea.
}
  \item{plot}{
Logical; should a diagnostic plot be made for human inspection? Top-left is the original \option{image}; top-right the median subtract sky polynomial; bottom-left the final sky subtract image (i.e. \option{image}[in] - \option{sky}[out]); bottom-right sky pixel density scaled by RMS (green dash is the Normal expectation, which is hopefully 1-1 in the negative domain).
}
  \item{\dots}{
Argument to pass to the 'Original' and 'Final' panels of the diagnostic plot when \option{plot} = TRUE.
}
}
\details{
Internally this uses the R \code{\link{lm}} and \code{\link{poly}} functions.

\code{profoundSkyPlane} is just a convenience function to setting \code{profoundSkyPoly} \option{degree} = 1 (i.e. a linear plane), since sometimes this confuses people.
}
\value{
\item{sky}{Numeric matrix; the predicted sky map that should be subtract from the user supplied \option{sky}.}
\item{skyRMS}{Numeric scalar; the calculated sky RMS. This is only indicative if the input \option{image} is from a single detector, in more complex cases a full map as provided by \code{\link{profoundMakeSkyGrid}} might be required.}
\item{lm_out}{List; the output of the \code{\link{lm}} run internally. This might be useful for additional analysis.}
\item{good_pix}{Integer matrix; mask of the 'good' pix, i.e. those valued 1 were used for analysis, and those valued 0 were not.}
\item{skyChiSq}{Integer scalar; the Reduced Chi-Square of the sky pixel distribution (the nearer to 1 the better the sky estimation).}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{lm}}, \code{\link{profoundMakeSkyGrid}}
}
\examples{
\dontrun{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat

profoundSkyPoly(image, degree=0, quancut=0.5, plot=TRUE) #under cut
profoundSkyPoly(image, degree=0, quancut=0.9, plot=TRUE) #about right
profoundSkyPoly(image, degree=0, quancut=0.999, plot=TRUE) #over cut

profoundSkyPoly(image, degree=1, quancut=0.5, plot=TRUE) #under cut
profoundSkyPoly(image, degree=1, quancut=0.9, plot=TRUE) #about right
profoundSkyPoly(image, degree=1, quancut=0.999, plot=TRUE) #over cut

profoundSkyPoly(image, degree=2, quancut=0.5, plot=TRUE) #under cut
profoundSkyPoly(image, degree=2, quancut=0.9, plot=TRUE) #about right
profoundSkyPoly(image, degree=2, quancut=0.999, plot=TRUE) #over cut

profoundSkyPoly(image, degree=3, quancut=0.5, plot=TRUE) #under cut
profoundSkyPoly(image, degree=3, quancut=0.9, plot=TRUE) #about right
profoundSkyPoly(image, degree=3, quancut=0.999, plot=TRUE) #over cut
}
}

\concept{ sky }
