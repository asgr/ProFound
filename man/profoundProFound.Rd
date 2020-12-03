\name{ProFound}
\alias{profoundProFound}
\alias{profound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound Source Detection
}
\description{
This is the highest level source detection function provided in \code{ProFit}, calculating both the initial segmentation map and reasonable estimates for the total flux apertures for each source in an automatic manner.
}
\usage{
profoundProFound(image = NULL, segim = NULL, objects = NULL, mask = NULL,  skycut = 1,
pixcut = 3, tolerance = 4, ext = 2, reltol = 0, cliptol = Inf, sigma = 1, smooth = TRUE,
SBlim, SBdilate = NULL, SBN100 = 100, size = 5, shape = "disc", iters = 6,
threshold = 1.05, magzero = 0, gain = NULL, pixscale = 1, sky = NULL, skyRMS = NULL,
redosegim = FALSE, redosky = TRUE, redoskysize = 21, box = c(100,100), grid = box,
skygrid_type = 'new', type = "bicubic", skytype = "median", skyRMStype = "quanlo",
roughpedestal = FALSE, sigmasel = 1, skypixmin = prod(box)/2, boxadd = box/2,
boxiters = 0, conviters = 100, iterskyloc = TRUE, deblend = FALSE, df = 3, radtrunc = 2,
iterative = FALSE, doclip = TRUE, shiftloc = FALSE, paddim = TRUE, header,
verbose = FALSE, plot = FALSE, stats = TRUE, rotstats = FALSE, boundstats = FALSE,
nearstats = boundstats, groupstats = boundstats, group = NULL, groupby = 'segim_orig',
offset = 1, haralickstats = FALSE, sortcol = "segID", decreasing = FALSE,
lowmemory = FALSE, keepim = TRUE, watershed = 'ProFound', pixelcov = FALSE,
deblendtype = 'fit', psf = NULL, fluxweight = 'sum', convtype = 'brute',
convmode = 'extended', fluxtype = 'Raw', app_diam = 1, Ndeblendlim = Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. If \option{image} is a list as created by \code{readFITS}, \code{read.fits} of \code{\link{magcutoutWCS}} then the image part of these lists is passed to \option{image} and the correct header part is passed to \option{header}. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This matrix *must* be the same dimensions as \option{image} if supplied. If this option is used then \code{profoundProFound} will not compute its initial segmentation map using \code{\link{profoundMakeSegim}}, which is then dilated. Instead it will use the one passed through \option{segim}.
}
  \item{objects}{
Boolean matrix (1,0); optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}. If provided then this is used to initially mask pixels for determining the correct \option{sky} when \option{sky}=NULL.
}
  \item{mask}{
Boolean matrix or integer scalar (1,0); optional, parts of the image to mask out (i.e. ignore). If a matrix is provided, this matrix *must* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis. if a scalar is provided it indicates the exact \option{image} values that should be treated as masked (e.g. by setting masked pixels to 0 or -999). The latter achieves the same effect as setting masked \option{image} pixels to NA, but allows for the fact not all programs can produce \code{R} legal NA values.
}
  \item{skycut}{
Numeric scalar; the lowest threshold to make on the \option{image} in units of the skyRMS. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbours, which is the highest. The range 1-5 offers decent results usually. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{ext}{
Numeric scalar; radius of the neighbourhood in pixels for the detection of neighbouring objects. Higher value smooths out small objects. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{reltol}{
Numeric scalar; only relevant for \option{watershed}='ProFound'. A modifier to the \option{tolerance}, modifying it by the ratio of the segment peak flux divided by the saddle point flux to the power \option{reltol}. The default means the \option{reltol} has no effect since this modifier becomes 1. A larger value of \option{reltol} means segments are more aggressively merged together. Can be (and often should be in practice) negative. The effect of using \option{reltol} and setting to negative is that the central brighter parts of galaxies are kept together in a single segment, and deblending is more common on the outskirts (where it should have less effect on the overall flux). The principle is that we need to be very confident a bright source needs to be split near its peak flux, but can be more aggressive in the outskirts.
}
  \item{cliptol}{
Numeric scalar; only relevant for \option{watershed}='ProFound'. If (\option{image}-\option{sky})/\option{skyRMS} is above this level where segments touch then they are always merged, regardless of other criteria. When thinking in terms of sky RMS, values between 20-100 are probably appropriate for merging very bright parts of stars back together in optical data.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE. Passed to \code{\link{profoundMakeSegim}}.
}
  \item{smooth}{
Logical; should smoothing be done on the target \option{image}? Passed to \code{\link{profoundMakeSegim}}. If present, this will use the \code{imblur} function from the \code{imager} package. Otherwise it will use the \code{gblur} function from the \code{EBImage} package with a warning. These functions are very similar in output, but not strictly identical.
}
  \item{SBlim}{
Numeric scalar; the mag/asec^2 surface brightness threshold to apply. This is always used in conjunction with \option{skycut}, so set \option{skycut} to be very large (e.g. Inf) if you want a pure surface brightness threshold for the segmentation. \option{magzero} and \option{pixscale} must also be present for this to be used. If a matrix is provided, this matrix *must* be the same dimensions as \option{image}. Passed to \code{\link{profoundMakeSegim}}. If set to 'get' then the output \option{SBlim} is calculated from the \option{skyRMS} and \option{magzero} automatically.
}
  \item{SBdilate}{
Numeric scalar; how many surface brightness mags beyond the sky RMS to push the dilation process. Default is NULL, meaning this is not triggered. The extra dilation logic is that if the new mean surface brightness in the dilated annulus is above the sky surface brightness + \option{SBdilate} and we are adding more than \option{SBN100} pixels then we will continue to dilate. This works well to capture extra flux in very large extended structures, where the fractional flux might not be changing much, but where we can safely keep dilating. Sensible values are between 0-2 (i.e. with 100 pixels we might be able to reliably push 2 mags fainter than the nominal sky surface brightness).
}
  \item{SBN100}{
Integer scalar; the number of new annulus pixels in our dilated segment required to trigger the \option{SBdilate} criteria.
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}. Passed to \code{\link{profoundMakeSegimDilate}}.
}
  \item{shape}{
Character scalar; the shape of the dilation kernel. See \code{makeBrush}. Passed to \code{\link{profoundMakeSegimDilate}}.
}
  \item{iters}{
Integer scalar; the maximum number of curve of growth dilations that should be made. This needs to be large enough to capture all the flux for sources of interest, but increasing this will increase the computation time for \code{profoundProFound}.  If this is set to 0 then the undilated \option{segim} image, whether provided or computed internally via \code{\link{profoundMakeSegim}}, will be used instead.
}
  \item{threshold}{
Numeric scalar; After the curve of growth dilations, \option{threshold} is the relative change of the converging property (see \option{converge}) that flags convergence. If consecutive iterations have a relative difference within this ratio then the dilation is stopped, and this iteration is used to define the segmentation of the object. The effect of this is that different objects will be dilated for a different number of iterations. Usually fainter sources require more. Note from v1.4: covergence has to be monotonically decreasing as dilation occurs, i.e. the flux can grow from x1.4 to x1.1 (next iteration), but not x1.4 to x2. Dilations with increasing growth rate are ignored.
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
User provided estimate of the absolute sky level. If this is not provided then it will be computed internally using \code{\link{profoundMakeSkyGrid}}. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{skyRMS}{
User provided estimate of the RMS of the sky. If this is not provided then it will be computed internally using \code{\link{profoundMakeSkyGrid}}. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).
}
  \item{redosegim}{
Logical; should the segmentation map be modified based using the interim "better sky"?  This means pixels falling below the new \option{skycut} would be excluded from the final segmentation map. This is usually only required if the sky subtraction was radically poor and complex in the first place. Will be forced to FALSE if the user supplies a segmentation map. If the user wants to flag object pixels they should pass it to \option{objects}.
}
  \item{redosky}{
Logical; should the sky and sky RMS grids be re-computed using the final segmentation map? This uses \code{\link{profoundMakeSkyGrid}} to compute the sky and sky RMS grids. If \option{redosky}=TRUE then the output will include the aggressively masked \option{objects_redo} image, if \option{redosky}=FALSE then \option{objects_redo} will be NA.
}
  \item{redoskysize}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels to apply to the \option{object} mask before performing the initial and final aggressively masked sky estimates (the latter is only relevant if \option{redosky}=TRUE). Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}. Dilation is done by \code{\link{profoundMakeSegimDilate}}. If \option{redosky}=TRUE, the final dilated \option{objects} mask is returned as \option{objects_redo}. As a rule of thumb you probably want ~50\% of your image pixels to be masked as objects, much more than this and you might not be able to sample enough sky pixels, much more less and the sky estimates might be biased by object flux in the wings.
}
  \item{box}{
Integer vector; the dimensions of the box car filter to estimate the sky with. For convenience, if length 1 then both dimensions of \option{box} used internally are assumed to equal the specified \option{box}. I.e. 200 would be interpreted as c(200,200). Dependent default arguments (\option{grid}, \option{boxadd} and \option{skypixmin}) are updated sensibly.
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
  \item{roughpedestal}{
Logical; when the initial "rough sky" is computed, should only a pedestal (based on the median of the sky map) be used for the sky? This is a good option if the image is known to contain a *very* large (many times the \option{box} size) galaxy that might otherwise be over subtracted by the initial rough sky map.
}
  \item{sigmasel}{
Numeric scalar; the quantile to use when trying to estimate the true standard-deviation of the sky distribution. If contamination is low then the default of 1 is about optimal in terms of S/N, but you might need to make the value lower when contamination is very high.
}
  \item{skypixmin}{
Numeric scalar; the minimum number of sky pixels desired in our cutout. The default is that we need half the original number of pixels in the \option{box} to be sky.
}
  \item{boxadd}{
Integer vector; the dimensions to add to the \option{box} to capture more pixels if \option{skypixmin} has not been achieved. By default this is set to be the same as the \option{box}/2.
}
  \item{boxiters}{
Integer scalar; the number of \option{box}+\option{boxadd} iterations to attempt in order to capture \option{skypixmin} sky pixels. The default means the box will not be grown at all.
}
  \item{conviters}{
Integer scalar; number of iterative sky convergence steps when \option{skytype} = 'converge' and/or \option{skyRMStype} = 'converge'. 
}
  \item{iterskyloc}{
Logical; should the last segment dilation be used to estimate a local sky value? If this is TRUE then this estimate is used to determine flux convergence for the dilations. It will also be used to return the additional \option{skyseg_mean} column in the segstats output provided. If TRUE then one additional dilation is made compared to the specified \option{iters}.
}
  \item{deblend}{
Logical; should segment flux be deblended using \code{\link{profoundFluxDeblend}} and these columns appended to the end of the output \option{segstats}?  
}
  \item{df}{
Integer scalar; degrees of freedom for the non-parametric spline fitting. Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{radtrunc}{
Numeric scalar; the maximum allowed radius beyond the edge-most segment pixel to consider when deblending. Keeping this low (1-3) ensures segments do not gather flux from very distant regions of the group. Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{iterative}{
Logical; should each segment profile fit be subtracted as it goes along? TRUE tends to remove the pedestal from a large galaxy that has faint objects embedded on top. Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{doclip}{
Logical; should the unmasked non-object pixels used to estimate to local sky value be further sigma-clipped using \code{\link{magclip}}? Whether this is used or not is a product of the quality of the objects extraction. If all detectable objects really have been found and the dilated objects mask leaves only apparent sky pixels then an advanced user might be confident enough to set this to FALSE. If in doubt, leave as TRUE.
}
	\item{shiftloc}{
Logical; should the cutout centre for the sky shift from \option{loc} of the desired \option{box} size extends beyond the edge of the image? (See \code{\link{magcutout}} for details).
}
	\item{paddim}{
Logical; should the cutout be padded with image data until it meets the desired \option{box} size (if \option{shiftloc} is true) or padded with NAs for data outside the image boundary otherwise? (See \code{\link{magcutout}} for details).
}
  \item{header}{
Full FITS header in table or vector format. If this is provided then the segmentations statistics table will gain \option{RAcen} and \option{Decen} coordinate outputs. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since a big \option{image} can take a long time to run, you might want to monitor progress.  
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{stats}{
Logical; should statistics on the segmented objects be returned using \code{\link{profoundSegimStats}}? If \option{magzero} and \option{pixscale} have been provided then some of the outputs are computed in terms of magnitude and mag/asec^2 rather than flux and flux/pix^2 (see Value).
}
  \item{rotstats}{
Logical; if TRUE then the \option{asymm}, \option{flux_reflect} and \option{mag_reflect} are computed, else they are set to NA. This is because they are very expensive to compute compared to other photometric properties.
}
  \item{boundstats}{
Logical; if TRUE then various pixel boundary statistics are computed (\option{Nedge}, \option{Nsky}, \option{Nobject}, \option{Nborder}, \option{edge_frac}, \option{edge_excess} and \option{FlagBorder}). If FALSE these return NA instead (saving computation time). Note by construction \option{Nedge} = \option{Nobject} + \option{Nsky} + \option{Nborder}. If you want to adjust specifically for \option{Nmask} then \option{Nsky} = \option{Nsky} - \option{Nmask}.
}
  \item{nearstats}{
Logical; if TRUE then the IDs of nearby segments is calculated via \code{\link{profoundSegimNear}} and output to the returned object \option{near}. By default this option is linked to \option{boundstats}, i.e. it is assumed if you want boundary statistics then you probably also want nearby object IDs returned.
}
  \item{groupstats}{
Logical; if TRUE then the IDs of grouped dilated segments (based on the output \option{segim}) is calculated via \code{\link{profoundSegimGroup}} and output to the list object \option{group}. By default this option is linked to \option{boundstats}, i.e. it is assumed if you want boundary statistics then you probably also want grouped object information returned. If \option{stats}=TRUE is also set then this flag will also create the \option{groupstats} output of photometric properties of the groups.
}
  \item{group}{
List; you can pass in the output from \code{\link{profoundSegimGroup}} directly, meaning groups will not be re-computed internally. This might be useful for speed in certain matched photometry applications.
}
  \item{groupby}{
Character scalar; How should the grouped segmentation map be formed that will be used to produce the \option{groupstats} output? Options are either via 'segim' or 'segim_orig'. 'segim' will create more groups, 'segim_orig' will have less.
}
  \item{offset}{
Integer scalar; the distance to offset when searching for nearby segments (used in both \code{\link{profoundSegimStats}} and \code{\link{profoundSegimNear}}).  
}
  \item{haralickstats}{
Logical; if TRUE then the Haralick texture statistics are computed using the \code{EBImage} function \code{computeFeatures.haralick}. For more detail see the original paper: http://haralick.org/journals/TexturalFeatures.pdf, and a useful online \code{EBImage} document: http://earlglynn.github.io/RNotes/package/EBImage/Haralick-Textural-Features.html.
}
  \item{sortcol}{
Character scalar; name of the output column that the returned segmentation statistics data.frame should be sorted by (the default is segID, i.e. segment order). See below for column names and contents.
}
  \item{decreasing}{
Logical; if FALSE (default) the segmentation statistics data.frame will be sorted in increasing order, if TRUE the data.frame will be sorted in decreasing order.
}
\item{lowmemory}{
Logical; if TRUE then a low memory mode of ProFound will be used. This limits the large \option{image} pixel matched outputs to just \option{segim}, with \option{segim_orig}, \option{objects} and \option{objects_redo} set to NULL, and \option{sky} and \option{skyRMS} set to 0. Internally the sky and skyRMS are used as normal for flux estimates, but they are removed as soon as possible within the function in order to free up memory. 
}
\item{keepim}{
Logical; if TRUE then the input \option{image} and \option{mask} matrices are passed through to the image output of the function. If FALSE then this is set to NULL.
}
  \item{watershed}{
Character scalar; the funciton to use to achieve the watershed deblend. Allowed options are 'EBImage' for EBImage::watershed, and 'ProFound' for the new \code{Rcpp} implementation included with the \code{ProFound} package.
}
  \item{pixelcov}{
Logical, should pixel covariance be considered when computing errors? TRUE uses \code{\link{profoundPixelCorrelation}} and the \option{cor_err_func} output which is then passed into \code{\link{profoundSegimStats}}. You can also directly pass it the \option{cor_err_func} function that you have already computed for the type of data being used. This means you do not need to recompute the correlation function for every \option{image} (which would be quite slow).
}
  \item{deblendtype}{
Character scalar, either 'fit' (default, where the segments in the image are approximately fitted for deblending) or 'psf' (where a PSF must be provided to argument \option{psf}). Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{psf}{
Numeric matrix; must be provided if \option{deblendtype}='psf'. Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{fluxweight}{
Numeric scalar; only relevent when \option{deblendtype}='psf'. Either 'sum' (where the sum of the current segment weights the deblend solution), 'peak' (where only the peak pixel flux in the segment it used) or 'none' (no additional weighting is used). Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{convtype}{
Character scalar, only relevent when \option{deblendtype}='psf'. Specifies the type on convolver to use. Available options are displayed by running \code{profitAvailableConvolvers}, but usually one of 'brute' or 'fftw'. The latter tends to be fastest when the supplied \option{psf} is much smaller than the \option{image} (factor greater than 4 in each dimension). As the sizes become comparable, fftw becomes faster since it scales better with the product of \option{psf} and \option{image} pixels.
}
  \item{convmode}{
Character scalar, only relevent when \option{deblendtype}='psf'. Either 'extended' (the whole segment is convolved with the target PSF when determining the deblend) or 'psf' (the provided PSF is centred at the peak flux for deblending). Only relevent when \option{deblendtype}='psf'. The former makes more sense for well resolved images (say optical) where the main issue in the deblending is the overlapping of resolved flux components. The latter works better in the regime where the image is barely resolved beyond the PSF, and this dominates the uncertainty of the deblend. 
}
  \item{fluxtype}{
Character scaler; specifies whether fluxes will be output in Jansky / MicroJansky ('Jansky' / 'microjansky'), or in raw/unscaled image ADUs ('Raw' / 'ADU' / 'ADUs', the default). You can only use Jansky / MicroJansky if the specified \option{magzero} gets the data into the AB system, else the fluxes will not be Jansky.
}
  \item{app_diam}{
Numeric scalar; the diameter in arc seconds to use for pseudo aperture photometry. This will use the appropriate pixel scale to convert the aperture into image units. The psuedo aperture photometry is output to columns \option{flux_app} and \option{mag_app} in \option{segstats}.
}
  \item{Ndeblendlim}{
Integer scalar; the limit for the number of pixels to consider in a deblending complex (Ngroup [number of segments in the group] x Npix [number of pixels in the group]). You might want to set this to a value similar to allowable machine memory (1e8 - 1e9 often) just to avoid extreme cases (e.g. large stars with lots of pixels and lots of segments). Only relevant if \option{deblend}=TRUE, see \code{\link{profoundFluxDeblend}}.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This high level function is both a source detection and a segmented aperture growing function. The latter is achieved through consecutive dilation and flux measurement operations. It is not super fast, but it is designed to be fairly robust and fast enough for most use cases.

\code{profoundProFound} initially makes a segmentation map using the \code{\link{profoundMakeSegim}} function. It then makes repeated dilations and flux measurements of this segmentation map using \code{\link{profoundMakeSegimDilate}}, and calculates the convergent flux segment for each source. These are combined to make a final segmentation map with associated source statistics (if requested).

The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut}, \option{pixcut}, \option{tolerance}, \option{sigma}, \option{ext}.

\code{\link{profoundMakeSegimDilate}} is similar in nature to the pixel growing \code{objmask} routine in \code{IRAF} (see the \option{ngrow} and \option{agrow} description at \url{http://iraf.noao.edu/projects/ccdmosaic/objmasks/objmasks.html}). This similarity was discovered after implementation, but it is worth noting that the higher level curve of growth function \code{profoundProFound} is not trivially replicated by other astronomy tools.
}
\value{
An object list of class 'profound' containing:

  \item{segim}{Integer matrix; the dilated and converged segmentation map matched pixel by pixel to \option{image}.}
  \item{segim_orig}{Integer matrix; the pre-dilated segmentation map matched pixel by pixel to \option{image}.}
  \item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{objects_redo}{Logical matrix; the dilated object map matched pixel by pixel to \option{image}. See \option{redosky} and \option{redoskysize}. Can be used as a mask in various other functions that require objects to be masked out.}
  \item{sky}{The estimated sky level of the \option{image}.}
  \item{skyRMS}{The estimated sky RMS of the \option{image}.}
  \item{image}{The input \option{image} matrix if \option{keepim}=TRUE, else NULL.}
  \item{mask}{The input or computed \option{mask} matrix if \option{keepim}=TRUE, else NULL.}
  \item{segstats}{If \option{stats}=TRUE then contains the output of \code{\link{profoundSegimStats}} run using the final \option{segim} (see below), otherwise NULL.}
  \item{Nseg}{The total number of segments extracted (dim(segstats)[1]).}
  \item{near}{If \option{nearstats}=TRUE then contains the output of \code{\link{profoundSegimNear}}.}
  \item{group}{If \option{groupstats}=TRUE then contains the output of \code{\link{profoundSegimGroup}}.}
  \item{groupstats}{If \option{groupstats}=TRUE and \option{stats}=TRUE then contains the output of \code{\link{profoundSegimStats}} run using the final \option{group$groupim} (see below), otherwise NULL.}
  \item{header}{The header provided, if missing this is NULL.}
  \item{SBlim}{The surface brightness limit of detected objects. Requires at least \option{magzero} to be provided and \option{skycut}>0, else NULL.}
  \item{magzero}{The assumed magnitude zero point. This is relevant to various outputs returned by the segmentation statistics.}
  \item{dim}{The dimensions of the processed image.}
  \item{pixscale}{The assumed pixel scale. This is relevant to various outputs returned by the segmentation statistics.}
  \item{gain}{The assumed image gain (if NULL it was not used). This is relevant to various outputs returned by the segmentation statistics.}
  \item{imarea}{The area of the provided image in degrees squared.}
  \item{skyLL}{The log-likelihood of negative sky pixel Chi-Sq distribution (the larger the better the sky estimation).}
  \item{call}{The original function call.}
  \item{date}{The date, more specifically the output of \code{\link{date}}.}
  \item{time}{The elapsed run time in seconds.}
  \item{ProFound.version}{The version of \code{ProFound} run, more specifically the output of \code{\link{packageVersion}}('ProFound').}
  \item{R.version}{The version of \code{R} run, more specifically the output of \code{\link{R.version}}.}
  
If \option{stats}=TRUE then the function \code{\link{profoundSegimStats}} is called and the \option{segstats} part of the returned list will contain a data.frame with columns (else NULL):

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
  \item{flux}{Total flux (calculated using \option{image}-\option{sky}) in ADUs or Jansky}
  \item{mag}{Total flux converted to mag using \option{magzero}}
  \item{flux_app}{Pseudo aperture (as specified by \option{app_rad}) flux (calculated using \option{image}-\option{sky}) in ADUs or Jansky}
  \item{mag_app}{Pseudo aperture (as specified by \option{app_rad}) flux converted to mag using \option{magzero}}
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
  \item{FPlim}{Approximate  false-positive significance limit below which one such source might appear spuriously on an image this large}
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
  \item{iter}{The iteration number when the source was flagged as having convergent flux}
  \item{origfrac}{The ratio between the final converged flux and the initial \code{\link{profoundMakeSegim}} iso-contour estimate}
  \item{Norig}{Number of pixels in the non-dilated (i.e. original) segment. This will be >= \option{pixcut} by construction.}
  \item{skyseg_mean}{This only has a value if \option{iterskyloc}=TRUE (otherwise NA). If provided it represents the mean sky for the outer dilation annulus created when dilating, i.e. a very local estimate of the sky. To use it the fluxes and mags provided need to have the current sky (\option{sky_sum}) added back on and the new sky (\option{skyseg_mean} x \option{N100}) subtracted.}
  \item{flag_keep}{A suggested flag for selecting good objects. Objects flagged FALSE have hit the iteration limit and have grown their flux by more than the median for all objects at the iteration limit. This tends to suggest a problem with the sky in the location of the segment.}
}
\author{
Aaron Robotham
}
\references{
Robotham A.S.G., et al., 2018, MNRAS, 476, 3137
Haralick R.M., et al., 1973, IEEE, SMC-3 (6), 610
}
\seealso{
\code{\link{profoundMakeSegim}}, \code{\link{profoundMakeSegimDilate}}, \code{\link{profoundMakeSegimExpand}}, \code{\link{profoundMakeSegimPropagate}}, \code{\link{profoundSegimStats}}, \code{\link{profoundSegimPlot}}, \code{\link{profoundMultiBand}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, magzero=30, verbose=TRUE, plot=TRUE)

# You can check to see if the final objects mask is aggressive enough. Notice the halos
# surrounding bright sources when just using the objects mask.

temp=image$imDat
temp[profound$objects>0]=0
magimage(temp)
temp=image$imDat
temp[profound$objects_redo>0]=0
magimage(temp)

magplot(profound$segstats[,c("R50","SB_N90")], log='x', grid=TRUE)
magplot(profound$segstats[,c("R50","SB_N90")], log='x', grid=TRUE)

magplot(profound$segstats[,c("flux","origfrac")], log='x', grid=TRUE)

## An example of a large galaxy:

VST_r=readFITS(system.file("extdata", 'VST_r.fits', package="magicaxis"))

# Running on defaults results in the central galaxy subtracting itself:
plot(profoundProFound(VST_r))

# Setting boxters=2 fixes things nicely:
plot(profoundProFound(VST_r, boxiters=2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ Segmentation }% use one of  RShowDoc("KEYWORDS")
\concept{ Detection }% __ONLY ONE__ keyword per line
