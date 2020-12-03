\name{profoundMakeSegim}
\alias{profoundMakeSegim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Watershed Image Segmentation
}
\description{
A mid level utility to achieve decent quality image segmentation. It uses a mixture of image smoothing and watershed segmentation propagation to identify distinct objects for use in, e.g., \code{profitSetupData} (where the \option{segim} list item output of \code{profoundMakeSegim} would be passed to the \option{segim} input of \code{profitSetupData}). For most user cases, people should use the higher level \code{\link{profoundProFound}} function.
}
\usage{
profoundMakeSegim(image = NULL, mask = NULL, objects = NULL, skycut = 1, pixcut = 3,
tolerance = 4, ext = 2, reltol = 0, cliptol = Inf, sigma = 1, smooth = TRUE, SBlim = NULL,
magzero = 0, gain = NULL, pixscale = 1, sky = NULL, skyRMS = NULL, header = NULL,
verbose = FALSE, plot = FALSE, stats = TRUE, rotstats = FALSE, boundstats = FALSE,
offset = 1, sortcol = "segID", decreasing = FALSE, watershed = 'ProFound', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{skycut}{
Numeric scalar; the lowest threshold to make on the \option{image} in units of the sky RMS. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of sky RMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbours, which is the highest. The range 1-5 offers decent results usually. This is passed to \option{tolerance} in 'EBImage', or \option{abstol} in 'ProFound' (see \option{watershed}).
}
  \item{ext}{
Numeric scalar; radius of the neighbourhood in pixels for the detection of neighbouring objects. Higher value smooths out small objects.
}
  \item{reltol}{
Numeric scalar; only relevant for \option{watershed}='ProFound'. A modifier to the \option{tolerance}, modifying it by the ratio of the segment peak flux divided by the saddle point flux to the power \option{reltol}. The default means the \option{reltol} has no effect since this modifier becomes 1. A larger value of \option{reltol} means segments are more aggressively merged together.
}
  \item{cliptol}{
Numeric scalar; only relevant for \option{watershed}='ProFound'. If (\option{image}-\option{sky})/option{skyRMS} is above this level where segments touch then they are always merged, regardless of other criteria. When thinking in terms of sky RMS, values between 20-100 are probably appropriate for merging very bright parts of stars back together in optical data.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}? If present, this will use the \code{imblur} function from the \code{imager} package. Otherwise it will use the \code{gblur} function from the \code{EBImage} package with a warning. These functions are very similar in output, but not strictly identical.
}
  \item{SBlim}{
Numeric scalar; the mag/asec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used.
}
  \item{magzero}{
  Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{gain}{
Numeric scalar; the gain (in photo-electrons per ADU). This is only used to compute object shot-noise component of the flux error (else this is set to 0).
}
  \item{pixscale}{
  Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output is in terms of pixels, otherwise it is in arcseconds. If provided along with \option{magzero} then the flux and surface brightness outputs will represent magnitudes and mag/asec^2.
}
  \item{sky}{
User provided estimate of the absolute sky level. If this is not provided then it will be computed internally using \code{\link{profoundSkyEst}}. Can be a scalar (value uniformly applied) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed internally using \code{\link{profoundSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since big image can take a long time to run, you might want to monitor progress.  
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned? If \option{magzero} and \option{pixscale} have been provided then some of the outputs are computed in terms of magnitude and mag/asec^2 rather than flux and flux/pix^2 (see Value).
}
  \item{rotstats}{
Logical; if TRUE then the \option{asymm}, \option{flux_reflect} and \option{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
}
  \item{boundstats}{
Logical; if TRUE then various pixel boundary statistics are computed (\option{Nedge}, \option{Nsky}, \option{Nobject}, \option{Nborder}, \option{edge_frac}, \option{edge_excess} and \option{FlagBorder}). If FALSE these return NA instead (saving computation time).
}
  \item{offset}{
Integer scalar; the distance to offset when searching for nearby segments (used in \code{\link{profoundSegimStats}}).
}
  \item{sortcol}{
Character; name of the output column that the returned segmentation statistics data.frame should be sorted by (the default is segID, i.e. segment order). See below for column names and contents.
}
  \item{decreasing}{
Logical; if FALSE (default) the segmentation statistics data.frame will be sorted in increasing order, if TRUE the data.frame will be sorted in decreasing order.
}
  \item{watershed}{
Character; the funciton to use to achieve the watershed deblend. Allowed options are 'EBImage' for EBImage::watershed, and 'ProFound' for the new \code{Rcpp} implementation included with the \code{ProFound} package.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
To use this function you will need to have \code{EBImage} installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

> source("http://bioconductor.org/biocLite.R")

> biocLite("EBImage")

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants).

The \code{profoundMakeSegim} function offers a high level internal to R interface for making quick segmentation maps. The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut}, \option{pixcut}, \option{tolerance}, \option{sigma}, \option{ext}.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{sky}{The estimated sky level of the \option{image}.}
  \item{skyRMS}{The estimated sky RMS of the \option{image}.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{header}{The header provided, if missing this is NULL.}
  \item{call}{The original function call.}
  
If \option{stats}=TRUE then the function \code{\link{profoundSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.frame else NULL (see \code{\link{profoundSegimStats}}).
}
\references{
See ?EBImage::watershed
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundMakeSegimExpand}}, \code{\link{profoundProFound}}, \code{\link{profoundSegimStats}}, \code{\link{profoundSegimPlot}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
segim=profoundMakeSegim(image, plot=TRUE)

#Providing a mask entirely removes regions of the image for segmentation:
mask=matrix(0,dim(image)[1],dim(image)[2])
mask[1:80,]=1
profoundMakeSegim(image, mask=mask, plot=TRUE)

#Providing a previously created object map can sometimes help with detection (not here):
profoundMakeSegim(image, mask=mask, object=segim$objects, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segmentation }% use one of  RShowDoc("KEYWORDS")
\concept{ watershed }% __ONLY ONE__ keyword per line
