\name{profoundMakeSky}
\alias{profoundMakeSkyMap}
\alias{profoundMakeSkyGrid}
\alias{profoundMakeSkyBlur}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Sky Maps
}
\description{
The high level function computes the absolute sky and sky RMS level over an image at a scale defined locally by the \option{box} parameter. This coarse map can then be used to compute sky/skyRMS values for the local sky anywhere on an image. This function uses \code{\link{profoundSkyEstLoc}} to calculate the sky statistics for the subset boxcar regions.
}
\usage{
profoundMakeSkyMap(image = NULL, objects = NULL, mask = NULL, sky = 0,
  box = c(100,100), grid = box, skytype = "median", skyRMStype = "quanlo", sigmasel = 1,
  skypixmin = prod(box)/2, boxadd = box/2, boxiters = 0, conviters = 100, doclip = TRUE,
  shiftloc = FALSE, paddim = TRUE, cores = 1)
  
profoundMakeSkyGrid(image = NULL, objects = NULL, mask = NULL, sky = 0, 
  box = c(100,100), grid = box, skygrid_type = 'new', type = 'bicubic',
  skytype = "median", skyRMStype = "quanlo", sigmasel = 1, skypixmin = prod(box)/2,
  boxadd = box/2, boxiters = 0, conviters = 100, doclip = TRUE, shiftloc = FALSE,
  paddim = TRUE, cores = 1)

profoundMakeSkyBlur(image = NULL, objects = NULL, box = 100,
  sigma = mean(box)*(4/pi)/sqrt(12))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar (value uniformly applied) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted! The sky statistics will be computed with \option{sky} removed, and then the input \option{sky} is added back on. By this route it is possible to refine a best effort sky measurement iteratively. For images with extreme gradients this might achieve a better solution, and in particular a more accurate \option{skyRMS} map, where it might be systematically high if the gradient is strong across the scale of the requested \option{box} (since this will appear to be enhanced variance).
}
  \item{box}{
Integer vector; the dimensions of the box car filter to estimate the sky with.
}
  \item{grid}{
Integer vector; the resolution of the background grid to estimate the sky with. By default this is set to be the same as the \option{box}.
}
  \item{skygrid_type}{
Character scalar; either 'new' (the new ADACS C++ sky grid code) or 'old' (the older R based code as used for ProFound <= v1.10).  
}
  \item{type}{
Character scalar; either "bilinear" for bilinear interpolation or "bicubic" for bicubic interpolation. The former creates sharper edges, the later is guaranteed to be first order differentiable. As of \code{ProFound} v1.13.0 we use the IntpAkimaUniform2 version available from www.geometrictools.com under a Boost 1.0 license (replacing the older \code{Akima} package and ADACS implementations, the former being memory intensive and the latter caused some small numerical artefacts).
}
  \item{skytype}{
Character scalar; the type of sky level estimator used. Allowed options are 'median' (the default), 'mean', 'mode' and 'converge' (see \code{\link{profoundSkyEstLoc}} for an explanation of what these estimators do). In all cases this is the estimator applied to unmasked and non-object pixels. If \option{doclip}=TRUE then the pixels will be dynamically sigma clipped before the estimator is run. 
}
  \item{skyRMStype}{
Character scalar; the type of sky level estimator used. Allowed options are 'quanlo' (the default), 'quanhi', 'quanboth', 'sd' and 'converge' (see \code{\link{profoundSkyEstLoc}} for an explanation of what these estimators do). In all cases this is the estimator applied to unmasked and non-object pixels. If \option{doclip}=TRUE then the pixels will be dynamically sigma clipped before the estimator is run.
}
  \item{sigmasel}{
Numeric scalar; the quantile to use when trying to estimate the true standard-deviation of the sky distribution. If contamination is low then the default of 1 is about optimal in terms of S/N, but you might need to make the value lower when contamination is very high.
}
  \item{skypixmin}{
Numeric scalar; the minimum number of sky pixels desired in our cutout. The default is that we need half the original number of pixels in the \option{box} to be sky.
}
  \item{boxadd}{
Integer vector; the dimensions to add to the \option{box} to capture more pixels if \option{skypixmin} has not been achieved.
}
  \item{boxiters}{
Integer scalar; the number of \option{box}+\option{boxadd} iterations to attempt in order to capture \option{skypixmin} sky pixels. The default means the box will not be grown at all.
}
  \item{conviters}{
Integer scalar; number of iterative sky convergence steps when \option{skytype} = 'converge' and/or \option{skyRMStype} = 'converge'. 
}
  \item{doclip}{
Logical; should the unmasked non-object pixels used to estimate to local sky value be further sigma-clipped using \code{\link{magclip}}? Whether this is used or not is a product of the quality of the objects extraction. If all detectable objects really have been found and the dilated objects mask leaves only apparent sky pixels then an advanced user might be confident enough to set this to FALSE. If an doubt, leave as TRUE.
}
	\item{shiftloc}{
Logical; should the cutout centre for the sky shift from \option{loc} if the desired \option{box} size extends beyond the edge of the image? (See \code{\link{magcutout}} for details).
}
	\item{paddim}{
Logical; should the cutout be padded with image data until it meets the desired \option{box} size (if \option{shiftloc} is true) or padded with NAs for data outside the image boundary otherwise? (See \code{\link{magcutout}} for details).
}
  \item{cores}{
Integer scalar; how many cores should be used to calculate sky properties of the image. Given the overhead for parallel computing, this should probably only be above 1 for larger images. 
}
  \item{sigma}{
Numeric scalar; the standard deviation of the blur (positive) used for smoothing the sky pixels to make a sky map (see \code{isoblur} in \code{imager} package).
}
}
\details{
The matrix generated will have many fewer pixels than the original \option{image}, so it will need to be interpolated back onto the full grid by some mechanism in order to have 1-1 values for the sky and sky RMS.
}
\value{
\code{profoundMakeSkyMap} produces a list of two lists. The first (called \option{sky}) contains a list of x,y,z values for the absolute sky, and second (called \option{skyRMS}) contains a list of x,y,z values for the sky RMS. The grids returned are as coarse as the \option{grid} option provided.

\code{profoundMakeSkyGrid} produces a list of two lists. The first (called \option{sky}) is a matrix of values for the absolute sky. The second (called \option{skyRMS}) is a matrix of values for the absolute sky RMS. The image matrices returned are pixel matched to the input \option{image} using the specified interpolation scheme.

\code{profoundMakeSkyBlur} produces a matrix of the sky map (there is no sky RMS component). To work well sources must be well masked either with NA in the \option{image}, or using the \option{objects} matrix. Since this function can only provide the additive sky, it cannot be used in isolation to computer the necessary sky information (we need the sky RMS map for \code{\link{profoundProFound}}), but it might offer a better smoother sky (so replacing the sky component of \code{profoundMakeSkyGrid}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundSkyEst}}, \code{\link{profoundSkyEstLoc}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
magimage(image)
skymap = profoundMakeSkyMap(image, box=c(89,89))
magimage(skymap$sky)
magimage(skymap$skyRMS)

# Now again, masking out the known objects (will not help too much in this case):

segim=profoundMakeSegim(image, skycut=1.5, plot=TRUE)
segim_ex=profoundMakeSegimExpand(image, segim$segim, skycut=-Inf, plot=TRUE)

skymap=profoundMakeSkyMap(image, objects=segim_ex$objects, box=c(89,89))
magimage(skymap$sky, magmap=FALSE)
magimage(skymap$skyRMS, magmap=FALSE)

# We can bilinear interpolate this onto the full image grid:

skybil = profoundMakeSkyGrid(image, objects=segim_ex$objects, box=c(89,89),
type='bilinear')
magimage(skybil$sky, magmap=FALSE)
magimage(skybil$skyRMS, magmap=FALSE)

# Or we can bicubic interpolate this onto the full image grid:

skybic = profoundMakeSkyGrid(image, objects=segim_ex$objects, box=c(89,89), type='bicubic')
magimage(skybic$sky, magmap=FALSE)
magimage(skybic$skyRMS, magmap=FALSE)

# The differences tend to be at the edges:

magimage(skybil$sky-skybic$sky, magmap=FALSE)
magimage(skybil$skyRMS-skybic$skyRMS, magmap=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ sky }% use one of  RShowDoc("KEYWORDS")
\concept{ RMS }% __ONLY ONE__ keyword per line
