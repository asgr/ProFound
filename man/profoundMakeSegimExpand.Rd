\name{profoundMakeSegimExpand}
\alias{profoundMakeSegimExpand}
\alias{profoundMakeSegimDilate}
\alias{profoundDilate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation Map Expansion and Dilation
}
\description{
A high level utility to achieve decent quality image segmentation based on the expansion of a pre-existing segmentation map. It uses smoothing and local flux weighted comparisons to grow the current segmentation map so as to better identify distinct objects for use in, e.g., \code{profitSetupData}.
}
\usage{
profoundMakeSegimExpand(image = NULL, segim = NULL, mask = NULL, objects = NULL,
  skycut = 1, SBlim = NULL, magzero = 0, gain = NULL, pixscale = 1, sigma = 1,
  smooth = TRUE, expandsigma = 5, expand = "all", sky = NULL, skyRMS = NULL,
  keyvalues = NULL, verbose = FALSE, plot = FALSE, stats = TRUE, rotstats = FALSE,
  boundstats = FALSE, offset = 1, sortcol = "segID", decreasing = FALSE, ...)
  
profoundMakeSegimDilate(image = NULL, segim = NULL, mask = NULL, size = 9, shape = "disc",
  expand = "all", magzero = 0, gain = NULL, pixscale = 1, sky = 0, skyRMS = 0,
  keyvalues = NULL, verbose = FALSE, plot = FALSE, stats = TRUE, rotstats = FALSE,
  boundstats = FALSE, offset = 1, sortcol = "segID", decreasing = FALSE, nthreads = 1,
  ...)
  
profoundDilate(segim = NULL, size = 3, shape = 'disc', expand = "all", iters = 1,
  nthreads = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the image. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{skycut}{
Numeric scalar; the lowest threshold to make on the \option{image} in units of the skyRMS. Since we are restricted to expanding out pre-existing segmentation regions we can usually afford to make this value lower than the equivalent in \code{\link{profoundMakeSegim}}.
}
  \item{SBlim}{
Numeric scalar; the magnitude/arcsec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used.
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
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}? If present, this will use the \code{imblur} function from the \code{imager} package. Otherwise it will use the \code{gblur} function from the \code{EBImage} package with a warning. These functions are very similar in output, but not strictly identical.
}
  \item{expandsigma}{
Numeric scalar; standard deviation of the blur used when expanding out the \option{segim}. Roughly speaking if \option{skycut} is set to a low number (say -5) then the expansion will not be prevented by the local sky level and it will grow by the number of pixels specified by \option{expandsigma}.
}
  \item{expand}{
Integer vector; specifies which segmentation regions should be expanded by the segID integer reference. If left with the default \option{expand}='all' then all segments will be expanded.
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number.
}
  \item{shape}{
Character scalar; the shape of the dilation kernel. Options are 'box', 'disc', 'diamond', 'Gaussian', 'line'.
}
  \item{iters}{
Integer scalar; how many times the dilation kernel should be applied to the image.  
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{keyvalues}{
List; header values to be used for the WCS.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since big image can take a long time to run, you might want to monitor progress.  
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned?
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
  \item{nthreads}{
Integer scalar; the number of CPU threads to run on. Only valid if \code{ProFound} is built with OpenMP support.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
The basic behaviour of \code{profoundMakeSegimExpand} and  \code{profoundMakeSegimDilate} is to intelligently expand out image segments already identified by, e.g., \code{\link{profoundMakeSegim}}.

The \code{profoundMakeSegimExpand} defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut} (or \option{SBlim}), \option{expandsigma}, \option{sigma}.

\code{profoundMakeSegimDilate} is similar in nature to the pixel growing \code{objmask} routine in \code{IRAF} (see the \option{ngrow} and \option{agrow} description). This similarity was discovered after implementation, but it is worth noting that the higher level curve of growth function \code{\link{profoundProFound}} is not trivially replicated by other astronomy tools.

The main difference between \code{profoundMakeSegimExpand} and  \code{profoundMakeSegimDilate} is the former grows the expansion a bit more organically, whereas the latter always gives new pixels to the brighter object if in doubt. That said, \code{profoundMakeSegimDilate} often gives very similar solutions and runs about 10+ times faster, so might be the only option for larger images.

\code{profoundDilate} just does the simple task of dilating the target segmentation map with the target kernel defined by \option{size} and \option{shape}, with the returned object being an interger matrix with the same dimensions as the input \option{segim} (but with segments appropriately dilated).
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{sky}{The estimated sky level of the \option{image}. \code{profoundMakeSegimExpand} only).}
  \item{skyRMS}{The estimated sky RMS of the \option{image}. \code{profoundMakeSegimExpand} only).}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{keyvalues}{The keyvalues provided, if missing this is NULL.}
  \item{SBlim}{The surface brightness limit of detected objects. Requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL. \code{profoundMakeSegimExpand} only.}
  \item{call}{The original function call.}
  
If \option{stats}=TRUE then the function \code{\link{profoundSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.fram (else NULL).
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundMakeSegim}}, \code{\link{profoundProFound}}, \code{\link{profoundSegimStats}}, \code{\link{profoundSegimPlot}}
}
\examples{
\dontrun{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat
segim=profoundMakeSegim(image, plot=TRUE, skycut=2)
profoundMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=1)
profoundMakeSegimDilate(image, segim$segim, plot=TRUE)

#Some other examples:

profoundMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=0)
profoundMakeSegimExpand(image, segim$segim, plot=TRUE, skycut=-Inf, sigma=3)

profoundMakeSegimDilate(image, segim$segim, plot=TRUE, size = 15)
profoundMakeSegimDilate(image, segim$segim, plot=TRUE, size = 21)

#This expansion process is a *much* better idea then simply setting the original skycut
#to a low value like 1/0:
profoundMakeSegim(image, plot=TRUE, skycut = 1)
profoundMakeSegim(image, plot=TRUE, skycut = 0)

#Just the dilation
magimage(profoundDilate(segim$segim, size=11))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segmentation }% use one of  RShowDoc("KEYWORDS")
