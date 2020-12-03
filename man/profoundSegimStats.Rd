\name{profoundSegimInfo}
\alias{profoundSegimStats}
\alias{profoundSegimPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Segmentation Statistics
}
\description{
Basic summary statistics for image segments, e.g. aperture parameters, fluxes and surface brightness estimates. These might provide useful first guesses to ProFit fitting parameters (particularly \option{flux}, \option{axrat} and \option{ang}).
}
\usage{
profoundSegimStats(image = NULL, segim = NULL, mask = NULL, sky = NULL, skyRMS = NULL,
magzero = 0, gain = NULL, pixscale = 1, header = NULL, sortcol = "segID",
decreasing = FALSE, rotstats = FALSE, boundstats = FALSE, offset = 1, cor_err_func = NULL,
app_diam = 1)
profoundSegimPlot(image = NULL, segim = NULL, mask = NULL, sky = NULL, header = NULL,
col = rainbow(max(segim), end=2/3), profound = NULL, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
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
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{sortcol}{
Character; name of the output column that the returned segmentation statistics data.frame should be sorted by (the default is segID, i.e. segment order). See below for column names and contents.
}
  \item{decreasing}{
Logical; if FALSE (default) the segmentation statistics data.frame will be sorted in increasing order, if TRUE the data.frame will be sorted in decreasing order.
}
  \item{rotstats}{
Logical; if TRUE then the \option{asymm}, \option{flux_reflect} and \option{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
}
  \item{boundstats}{
Logical; if TRUE then various pixel boundary statistics are computed (\option{Nedge}, \option{Nsky}, \option{Nobject}, \option{Nborder} \option{Nmask}, \option{edge_frac}, \option{edge_excess} and \option{FlagBorder}). If FALSE these return NA instead (saving computation time). Note by construction \option{Nedge} = \option{Nobject} + \option{Nsky} + \option{Nborder}. If you want to adjust specifically for \option{Nmask} then \option{Nsky} = \option{Nsky} - \option{Nmask}.
}
  \item{offset}{
Integer scalar; the distance to offset when searching for nearby segments.  
}
  \item{col}{
Colour palette; the colours to map the segment IDs against. This is by default the magnitude using a rainbow palette, going from red for bright segments, via green, to blue for faint segments.
}
  \item{profound}{
List; object of class 'profound'. If this is provided then missing input arguments are taking directly from this structure. As an added convenience, you can assign the profound object directly to the \option{image} input.
}
  \item{cor_err_func}{
Function; the error function between \option{N100} (the number of pixels in the segment) and the relative flux error. Most likely the \option{cor_err_func} output of \code{\link{profoundPixelCorrelation}}.
}
  \item{app_diam}{
Numeric scalar; the diameter in arc seconds to use for pseudo aperture photometry. This will use the appropriate pixel scale to convert the aperture into image units. The psuedo aperture photometry is output to columns \option{flux_app} and \option{mag_app} in \option{segstats}.
}
  \item{add}{
Logical; should just the segment contours be added to the current image? This allows for complex colouring of different segments to be achieved by adding various overlays.  
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}.
}
}
\details{
\code{profoundSegimStats} provides summary statistics for the individual segments of the image, e.g. properties of the apertures, and the sum of the flux etc. This is used inside of \code{\link{profoundMakeSegim}} and \code{\link{profoundMakeSegimExpand}}, but it may be useful to use separately if manual modifications are made to the segmentation, or two segmentations (e.g. a hot and cold mode segmentation) need to be combined.

The interpretation of some of these outputs will depend a lot on the data being analysed, so it is for the user to decide on sensible next steps (e.g. using the outputs to select stars etc). One output of interest might be \option{flux_reflect}. This attempts to correct for missing flux where segments start colliding. This probably returns an upper limit to the flux since in some regions it can even be double counted if the two sources that have colliding segmentation maps are very close together and similar in brightness, so somewhere between \option{flux} and \option{flux_reflect} the truth probably lies. If you want a better estimate of the flux division then you should really be using the profiling routine of \code{ProFit}.

\code{profoundSegimPlot} is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot! If you provide a header or a list containing the iamge and header to \option{header} then it will be plotted with the WCS overlaid using \code{\link{magimageWCS}}, otherwise it will use \code{\link{magimage}}.
}
\value{
A data.frame with columns:

  \item{segID}{Segmentation ID, which can be matched against values in \option{segim}}
  \item{uniqueID}{Unique ID, which is fairly static and based on the xmax and ymax position}
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
  \item{flux_app}{Pseudo aperture (as specified by \option{Napp}) flux (calculated using \option{image}-\option{sky}) in ADUs or Jansky}
  \item{mag_app}{Pseudo aperture (as specified by \option{Napp}) flux converted to mag using \option{magzero}}
  \item{cenfrac}{Fraction of flux in the brightest pixel}
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
  \item{semimaj}{Weighted standard deviation along the major axis, i.e. the semi-major first moment, so ~2 times this would be a typical major axis Kron radius (always in units of pix)}
  \item{semimin}{Weighted standard deviation along the minor axis, i.e. the semi-minor first moment, so ~2 times this would be a typical minor axis Kron radius  (always in units of pix)}
  \item{axrat}{Axial ratio as given by min/maj}
  \item{ang}{Orientation of the semi-major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical)}
  \item{signif}{Approximate singificance of the detection using the Chi-Square distribution}
  \item{FPlim}{Approximate false-positive significance limit below which one such source might appear spuriously on an image this large}
  \item{flux_err}{Estimated total error in the flux for the segment}
  \item{mag_err}{Estimated total error in the magnitude for the segment}
  \item{flux_err_sky}{Sky subtraction component of the flux error}
  \item{flux_err_skyRMS}{Sky RMS component of the flux error}
  \item{flux_err_shot}{Object shot-noise component of the flux error (only if \option{gain} is provided)}
  \item{flux_err_cor}{Error component due to pixel correlation}
  \item{sky_mean}{Mean flux of the sky over all segment pixels}
  \item{sky_sum}{Total flux of the sky over all segment pixels}
  \item{skyRMS_mean}{Mean value of the sky RMS over all segment pixels}
  \item{Nedge}{Number of edge segment pixels that make up the outer edge of the segment}
  \item{Nsky}{Number of edge segment pixels that are touching sky}
  \item{Nobject}{Number of edge segment pixels that are touching another object segment}
  \item{Nborder}{Number of edge segment pixels that are touching the \option{image} border}
  \item{Nmask}{Number of edge segment pixels that are touching a masked pixel (note NAs in \option{image} are also treated as masked pixels)}
  \item{edge_frac}{Fraction of edge segment pixels that are touching the sky i.e. \option{Nsky}/\option{Nedge}, higher generally meaning more robust segmentation statistics}
  \item{edge_excess}{Ratio of the number of edge pixels to the expected number given the elliptical geometry measurements of the segment. If this is larger than 1 then it is a sign that the segment geometry is irregular, and is likely a flag for compromised photometry}
  \item{flag_border}{A binary flag telling the user which \option{image} borders the segment touches. The bottom of the \option{image} is flagged 1, left=2, top=4 and right=8. A summed combination of these flags indicate the segment is in a corner touching two borders: bottom-left=3, top-left=6, top-right=12, bottom-right=9.}
  
\code{profoundSegimPlot} is a simple function that overlays the image segments on the original \option{image}. This can be very slow for large numbers (1,000s) of segments because it uses the base \code{\link{contour}} function to draw the segments individually.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundProFound}}, \code{\link{profoundMakeSegim}}, \code{\link{profoundMakeSegimExpand}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, magzero=30, rotstats=TRUE)

print(profound$segstats)

#Note row 6 (the central galaxy) gains 0.05 mag of flux due to the missing flux when
#rotated through 180 degrees. The reflected value of 18.4 is closer to the full profile
#solution (~18.35) than the non-reflected flux (18.45).

profound$segim[35:55, 80:100]=max(profound$segim)+1
print(profoundSegimStats(image$imDat, segim=profound$segim, sky=profound$sky,
header=image$hdr))
profoundSegimPlot(image, profound$segim)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segmentation }% use one of  RShowDoc("KEYWORDS")
