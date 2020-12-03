\name{profoundSkyEstLoc}
\alias{profoundSkyEstLoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Sky in Subset of Pixels
}
\description{
Calculate the sky and sky RMS for a subset region of a larger image, as used in \code{\link{profoundMakeSkyMap}}.
}
\usage{
profoundSkyEstLoc(image = NULL, objects = NULL, mask = NULL, loc = dim(image)/2,
box = c(100, 100), skytype = "median", skyRMStype = "quanlo", sigmasel = 1,
skypixmin = prod(box)/2, boxadd = box/2, boxiters = 0, conviters = 100, doclip = TRUE,
shiftloc = FALSE, paddim = TRUE, plot = FALSE, ...)
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
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{loc}{
Integer vector; the [x,y] location where we want to estimate the sky and sky RMS.
}
  \item{box}{
Integer vector; the dimensions of the box car filter to estimate the sky with.
}
  \item{skytype}{
Character scalar; the type of sky level estimator used. Allowed options are 'median' (the default), 'mean', 'mode' and 'converge' (see Details for an explanation of what these estimators do). In all cases this is the estimator applied to unmasked and non-object pixels. If \option{doclip}=TRUE then the pixels will be dynamically sigma clipped before the estimator is run. 
}
  \item{skyRMStype}{
Character scalar; the type of sky level estimator used. Allowed options are 'quanlo' (the default), 'quanhi', 'quanboth', 'sd' and 'converge' (see Details for an explanation of what these estimators do). In all cases this is the estimator applied to unmasked and non-object pixels. If \option{doclip}=TRUE then the pixels will be dynamically sigma clipped before the estimator is run.
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
Logical; should the cutout center shift from \option{loc} if the desired \option{box} size extends beyond the edge of the image? (See \code{\link{magcutout}} for details).
}
	\item{paddim}{
Logical; should the cutout be padded with image data until it meets the desired \option{box} size (if \option{shiftloc} is true) or padded with NAs for data outside the image boundary otherwise? (See \code{\link{magcutout}} for details).
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This is a somewhat handy standalone utility function if you have a large image and want to check the quality and stability of the local sky and sky RMS.

Regarding \option{skytype}, the meaning of the 'median' and 'mean' options are obvious enough. The 'mode' is computed by running the data through \code{\link{density}} with the default options including automatuc selection of the appropriate smoothing band-width. The peak value of the smoothed density is then extracted, and the pixel value at this point is returned as the 'mode' sky estimator. The 'converge' sky uses a convergence scheme based on the estimated mean and variance of a truncated Normal distibution, where it attempts to maximise the likelihood of the population mode and standard deviation for the Normal sky.

Regarding \option{skyRMStype}, if you know that your contamination only comes from positive flux sources (e.g., astronomical data when trying to select sky pixels) then you should probably use the lower side to determine Normal statistics (quanlo). Similarly if the contamination is on the low side then you should use the higher side to determine Normal statistics (quanhi, but this is rare in astronomical data). If you believe the selected sky pixels to be unbiased then 'quanboth' uses both sides and will give you a more accurate estimator of the sky RMS. The 'sd' option is to use the standard-deviation, with the caveat that this is calculated around the esstimated sky level (of type specified by \option{skytype}) and not necessarily simply the mean (as it would be typically). The most common choices for \option{skyRMStype} will likely be 'quanlo' or 'sd'. The 'converge' sky uses a convergence scheme based on the estimated mean and variance of a truncated Normal distibution, where it attempts to maximise the likelihood of the population mode and standard deviation for the Normal sky.

There are many questions to think about when choosing the best combination of sky estimators. Have all detectable sources been robustly extracted and masked? Is the remaining contamintion due to background undetected sources or wing flux from foreground stars? The most significant choice to be made is whether to choose the more robust 'median' or the potentially biased 'mean'. The former makes sense if you think there might be detectable sources still contributing to your nominal sky pixels, the latter makes sense if the positive  flux of undetected sources is spread round the sky in an random but uniform manner. If you are very confident that your object mask represents all plausible sources then you might even want to set \option{doclip}=FALSE. The defaults behave in quite a safe manner and have resistance to unmasked objects being included in the sky pixels. Using different options (particularly \option{doclip}=FALSE and \option{skytype}) requires more advanced knowledge about the specific data being anlysed.

If the package \code{Rfast} is loaded, then this function will use the faster Rfast::med median function over the base median. This is about a factor 2-3 faster.
}
\value{
A length two vector where the first element is the sky and the second is the skyRMS.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSkyEst}}, \code{\link{profoundMakeSkyMap}}, \code{\link{profoundMakeSkyGrid}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
profoundSkyEstLoc(image, loc=c(20,20), box=c(40,40), plot=TRUE)$val
profoundSkyEstLoc(image, loc=c(40,20), box=c(40,40), plot=TRUE)$val
profoundSkyEstLoc(image, loc=c(60,20), box=c(40,40), plot=TRUE)$val
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ sky }% use one of  RShowDoc("KEYWORDS")
\concept{ RMS }% __ONLY ONE__ keyword per line
