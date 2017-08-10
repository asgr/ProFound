\name{profoundMakeSegim}
\alias{profoundMakeSegim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Watershed Image Segmentation
}
\description{
A high level utility to achieve decent quality image segmentation. It uses a mixture of image smoothing and watershed segmentation propagation to identify distinct objects for use in, e.g., \code{profitSetupData} (where the \option{segim} list item output of \code{profoundMakeSegim} would be parsed to the \option{segim} input of \code{profitSetupData}).
}
\usage{
profoundMakeSegim(image, mask, objects, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE,
pixcut = 5, skycut = 2, SBlim, magzero = 0, gain, pixscale = 1, sky, skyRMS, header,
verbose = FALSE, plot = FALSE, stats = TRUE, rotstats = FALSE,  boundstats = FALSE,
offset = 1, sortcol = "segID", decreasing = FALSE, ...)
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
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbours, which is the highest. The range 1-5 offers decent results usually.
}
  \item{ext}{
Numeric scalar; radius of the neighbourhood in pixels for the detection of neighbouring objects. Higher value smooths out small objects.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}?
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object.
}
  \item{skycut}{
Numeric scalar; the lowest threshold to make on the \option{image} in units of the skyRMS.
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
User provided estimate of the absolute sky level. If this is not provided then it will be computed internally using \code{\link{profoundSkyEst}}. Can be a scalar (value uniformly applied to full \option{sigma} map) or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
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
Character; name of the output column that the returned segmentation statistics data.frame should be sorted by. See below for column names and contents.
}
  \item{decreasing}{
Logical; if FALSE (default) the segmentation statistics data.frame will be sorted in increasing order, if TRUE the data.frame will be sorted in decreasing order.
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

The \code{profoundMakeSegim} function offers a high level internal to R interface for making quick segmentation maps. The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{tolerance}, \option{sigma}, \option{ext}, \option{skycut}, \option{pixcut}.
}
\value{
A list containing:

  \item{segim}{Integer matrix; the segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{sky}{The estimated sky level of the \option{image}.}
  \item{skyRMS}{The estimated sky RMS of the \option{image}.}
  \item{segstats}{If \option{stats}=TRUE this is a data.frame (see below), otherwise NULL.}
  \item{header}{The header provided, if missing this is NULL.}
  \item{SBlim}{The surface brightness limit of detected objects (requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL).}
  \item{call}{The original function call.}
  
If \option{stats}=TRUE then the function \code{\link{profoundSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.frame with columns (else NULL):

  \item{segID}{Segmentation ID, which can be matched against values in \option{segim}}
  \item{uniqueID}{Unique ID, which is fairly static and based on the xcen and ycen position}
  \item{xcen}{Flux weighted x centre}
  \item{ycen}{Flux weighted y centre}
  \item{xmax}{x position of maximum flux}
  \item{ymax}{y position of maximum flux}
  \item{RAcen}{Flux weighted degrees Right Ascension centre (only present if a \option{header} is provided)}
  \item{Deccen}{Flux weighted degrees Declination centre (only present if a \option{header} is provided)}
  \item{RAmax}{Right Ascension of maximum flux (only present if a \option{header} is provided)}
  \item{Decmax}{Declination of maximum flux (only present if a \option{header} is provided)}
  \item{sep}{Radial offset between the cen and max definition of the centre (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)}
  \item{flux}{Total flux (calculated using \option{image}-\option{sky}) in ADUs}
  \item{mag}{Total flux converted to mag using \option{magzero}}
  \item{N50}{Number of brightest pixels containing 50\% of the flux}
  \item{N90}{Number of brightest pixels containing 90\% of the flux}
  \item{N100}{Total number of pixels in this segment, i.e. contains 100\% of the flux}
  \item{R50}{Approximate elliptical semi-major axis containing 50\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)}
  \item{R90}{Approximate elliptical semi-major axis containing 90\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)}
  \item{R100}{Approximate elliptical semi-major axis containing 100\% of the flux (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)}
  \item{SB_N50}{Mean surface brightness containing brightest 50\% of the flux, calculated as \option{flux}*0.5/\option{N50} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2)}
  \item{SB_N90}{Mean surface brightness containing brightest 90\% of the flux, calculated as \option{flux}*0.9/\option{N90} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2)}
  \item{SB_N100}{Mean surface brightness containing all of the flux, calculated as \option{flux}/\option{N100} (if \option{pixscale} has been set correctly then this column will represent mag/asec^2. Otherwise it will be mag/pix^2)}
  \item{xsd}{Weighted standard deviation in x (always in units of pix)}
  \item{ysd}{Weighted standard deviation in y (always in units of pix)}
  \item{covxy}{Weighted covariance in xy (always in units of pix)}
  \item{corxy}{Weighted correlation in xy (always in units of pix)}
  \item{con}{Concentration, \option{R50}/\option{R90}}
  \item{asymm}{180 degree flux asymmetry (0-1, where 0 is perfect symmetry and 1 complete asymmetry)}
  \item{flux_reflect}{Flux corrected for asymmetry by doubling the contribution of flux for asymmetric pixels (defined as no matching segment pixel found when the segment is rotated through 180 degrees)}
  \item{mag_reflect}{\option{flux_reflect} converted to mag using \option{magzero}}
  \item{maj}{Weighted standard deviation along the major axis, i.e. the first moment, so ~2 times this would be a typical major axis Kron radius (always in units of pix)}
  \item{min}{Weighted standard deviation along the minor axis, i.e. the first moment, so ~2 times this would be a typical minor axis Kron radius  (always in units of pix)}
  \item{axrat}{Axial ratio as given by min/maj}
  \item{ang}{Orientation of the major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical)}
  \item{flux_err}{Estimated total error in the flux for the segment}
  \item{mag_err}{Estimated total error in the magnitude for the segment}
  \item{flux_err_sky}{Sky subtraction component of the flux error}
  \item{flux_err_skyRMS}{Sky RMS component of the flux error}
  \item{flux_err_shot}{Object shot-noise component of the flux error (only if \option{gain} is provided)}
  \item{sky_mean}{Mean flux of the sky over all segment pixels}
  \item{sky_sum}{Total flux of the sky over all segment pixels}
  \item{skyRMS_mean}{Mean value of the sky RMS over all segment pixels}
  \item{Nedge}{Number of edge segment pixels that make up the outer edge of the segment}
  \item{Nsky}{Number of edge segment pixels that are touching sky}
  \item{Nobject}{Number of edge segment pixels that are touching another object segment}
  \item{Nborder}{Number of edge segment pixels that are touching the \option{image} border}
  \item{Nmask}{Number of edge segment pixels that are touching a masked pixel (note NAs in \option{image} are also treated as masked pixels)}
  \item{edge_frac}{Fraction of edge segment pixels that are touching the sky i.e. \option{Nsky}\option{Nedge}, higher generally meaning more robust segmentation statistics}
  \item{edge_excess}{Ratio of the number of edge pixels to the expected number given the elliptical geometry measurements of the segment. If this is larger than 1 then it is a sign that the segment geometry is irregular, and is likely a flag for compromised photometry}
  \item{Flagborder}{A binary flag telling the user which \option{image} borders the segment touches. The bottom of the \option{image} is flagged 1, left=2, top=4 and right=8. A summed combination of these flags indicate the segment is in a corner touching two borders: bottom-left=3, top-left=6, top-right=12, bottom-right=9.}
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
\keyword{ segmentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ watershed }% __ONLY ONE__ keyword per line
