\name{profoundMultiBand}
\alias{profoundMultiBand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi Band ProFound Photometry
}
\description{
Run multiband \code{ProFound} photometry either with loaded data, or images on a local disk.
}
\usage{
profoundMultiBand(inputlist = NULL, dir = "", segim, mask, skycut = 1, pixcut = 3,
tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, iters_det = 6, iters_tot = 2,
detectbands = "r", multibands = c("u", "g", "r", "i", "z"), magzero = 0, gain = NULL,
bandappend = multibands, totappend = "t", colappend = "c", grpappend ='g', dotot = TRUE,
docol = TRUE, dogrp=TRUE, boundstats = TRUE, haralickstats = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputlist}{
A list of already loaded images. Typically of the type loaded in from FITS files by the \code{astro} package's \code{read.fits} function, or the \code{FITSio} package's \code{readFITS} function. If using the \option{inputlist} parameter the length of the list must be the same length as \option{multibands} (and the related parameters).
}
  \item{dir}{
If \option{inputlist} is left as NULL then \code{profoundMultiBand} will instead try to load in FITS images from the directory specified by \option{dir}. The images in the directory must have names like \option{multibands}[i].fits etc (so with the defaults names like u.fits and g.fits would be okay). Since \option{multibands} effectively specifies the file names much more complicated naming can be used and passed in, but it is also used by default for naming the catalogue column outputs, so shorter names/references are likely to be preferable there (i.e. mag_ut is simpler than mag_KiDS_VST_ut etc). This can be over-ridden by using \option{catappend}.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This matrix *must* be the same dimensions as the detection image/s if supplied. If this option is used then \code{profoundMultiBand} will not compute its initial segmentation map using \code{\link{profoundMakeSegim}}, which is then dilated. Instead it will use the one passed through \option{segim} and dilate this according to the \option{iters_det} argument (so set this to 0 if you want the \option{segim} to be used as is).
}
  \item{mask}{
Boolean matrix or integer scalar; optional, parts of the image to mask out (i.e. ignore). If a matrix is provided, this matrix *must* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis. if a scalar is provided it indicates the exact \option{image} values that should be treated as masked (e.g. by setting masked pixels to 0 or -999). The latter achieves the same effect as setting masked \option{image} pixels to NA, but allows for the fact not all programs can produce \code{R} legal NA values. 
}
  \item{skycut}{
Numeric scalar; the lowest threshold to make on the detection \option{image} in units of the skyRMS. Passed to \code{\link{profoundProFound}}.
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object. Passed to \code{\link{profoundProFound}}.
}
  \item{tolerance}{
Numeric scalar; the minimum height of the object in the units of skyRMS between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbours, which is the highest. The range 1-5 offers decent results usually. Passed to \code{\link{profoundProFound}}.
}
  \item{ext}{
Numeric scalar; radius of the neighbourhood in pixels for the detection of neighbouring objects. Higher value smooths out small objects. Passed to \code{\link{profoundProFound}}.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used when \option{smooth}=TRUE. Passed to \code{\link{profoundProFound}}.
}
  \item{smooth}{
Logical; should smoothing be done on the target detection \option{image}? Passed to \code{\link{profoundProFound}}. If present, this will use the \code{imblur} function from the \code{imager} package. Otherwise it will use the \code{gblur} function from the \code{EBImage} package with a warning. These functions are very similar in output, but not strictly identical.
}
  \item{iters_det}{
Integer scalar; the maximum number of curve of growth dilations that should be made to the detection image. This needs to be large enough to capture all the flux for sources of interest, but increasing this will increase the computation time for \code{profoundProFound}. If this is set to 0 then the undilated \option{segim} image, whether provided or computed internally via \code{\link{profoundMakeSegim}}, will be used instead.
}
  \item{iters_tot}{
Integer scalar; the maximum number of curve of additional growth dilations that should be made above the dilated detection segmentation map for multi band total colour photometry. This is only relevant if \option{dotot}=TRUE. This should not be set too high (and might even be 0) since the detection image should generally be fairly deep.
}
  \item{detectbands}{
Character vector; the names of the detection bands that will be stacked using \code{\link{profoundMakeStack}} and then analysed with the provided settings with \code{\link{profoundProFound}} to make a reference segmentation map for further multi band photometry. These bands must be present in \option{multibands}. Can be a scalar (i.e. a single band is used). If set to 'get' then it will use all legal FITS files in the target directory. If set to 'all' then it will use all \option{multibands} inputs.
}
  \item{multibands}{
Character vector; the names of the target multi band photometry images. If set to 'get' then it will use all legal FITS files in the target directory. If using the \option{inputlist} parameter the length of the list must be the same length as \option{multibands}. \option{magzero} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. If specified, \option{gain} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. If specified, \option{catappend} must either the same length and order as \option{multibands}.
}
  \item{magzero}{
Numeric vector; the magnitude zero point of the images being used. \option{magzero} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{gain}{
Numeric vector; the gain of the images being used. \option{gain} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{bandappend}{
Character vector; characters to be appended per band in the output multi band photometry catalogues. The default will create columns with names like mag_ut (total) and mag_uc (colour).
}
  \item{totappend}{
Character scalar; character to be appended in the output multi band total photometry catalogue (\code{cat_tot}). The default will create columns with names like mag_ut and R50_ut.
}
  \item{colappend}{
Character scalar; character to be appended in the output multi band colour photometry catalogue (\code{cat_col}). The default will create column with names like mag_uc and R50_uc.
}
  \item{grpappend}{
Character scalar; character to be appended in the grouped segment multi band total photometry catalogue (\code{cat_tot}). The default will create columns with names like mag_ug and R50_ug.
}
  \item{dotot}{
Logical; should dilated segment total photometry be computed for the bands specified in \option{multibands}? This will return closer to total magnitudes in all target bands.
}
  \item{docol}{
Logical; should non-dilated segment colour photometry be computed for the bands specified in \option{multibands}? This will return better colour magnitudes in all target bands (i.e. more accurate differences between bands) and will typically under-represent the total photometry.
}
  \item{dogrp}{
Logical; should group segment photometry be computed for the bands specified in \option{multibands}? This might be useful for re-assembling large galaxies that are broken up at a later date. \option{boundstats} must also be set to TRUE if \option{dogrp}=TRUE is set.
}
  \item{boundstats}{
Logical; if TRUE then various pixel boundary statistics are computed for the detection image and all multi-band photometry. See also \code{\link{profoundProFound}}. \option{boundstats} must also be set to TRUE if \option{dogrp}=TRUE is set.
}
  \item{haralickstats}{
Logical; if TRUE then the Haralick texture statistics are computed for the detection image. See also \code{\link{profoundProFound}}.
}
  \item{verbose}{
Logical; should verbose output from \code{\link{profoundProFound}} be displayed to the user? Basic information is always printed out, to help catch silly user errors quickly.
}
  \item{\dots}{
Further arguments to be passed to the stacked detection image run of \code{\link{profoundProFound}}.
}
}
\details{
This very high level function simplifies a sequence of function calls that we found users typically needed to make, but when scripted they were prone to mistakes and made multi band photometry scripts hard to maintain.

In the simplest sense this script runs \code{\link{profoundProFound}} on each detection band and uses this information to make a stacked image using \code{\link{profoundMakeStack}}. \code{\link{profoundProFound}} is then run on this stacked image to make a deep segmentation map. For good total photometry the segim object from this output is used, and allowed to further dilate to account for different observing conditions (i.e. PSFs). For good colour photometry the segim_orig object from this output is used. Only the \code{\link{profoundSegimStats}} output is kept for the target multi band images, so not all of the outputs from \code{\link{profoundProFound}} since this is usually unnecessary when operating in this mode, and creates a huge quantity of data.
}
\value{
An object list of class 'profoundmulti' containing:

\item{pro_detect}{The full output of \code{\link{profoundProFound}} for the detection image (of class 'profound').}
\item{cat_tot}{If \option{dotot}=TRUE, the dilated total photometry for the target bands. Effectively the output of \code{\link{profoundSegimStats} run on pro_detect$segim}.}
\item{cat_col}{If \option{docol}=TRUE, the non-dilated colour photometry for the target bands. Effectively the output of \code{\link{profoundSegimStats} run on pro_detect$segim_orig}.}
\item{cat_grp}{If \option{dogrp}=TRUE, the group segment photometry for the target bands. Effectively the output of \code{\link{profoundSegimStats} run on pro_detect$group$groupim}.}
\item{detectbands}{Character vector; the names of the detection bands used.}
\item{multibands}{Character vector; the names of the target multi band photometry images used.}
\item{call}{The original function call.}
\item{date}{The date, more specifically the output of \code{\link{date}}.}
\item{time}{The elapsed run time in seconds.}
}
\references{
Robotham A.S.G., et al., 2018, MNRAS, 476, 3137
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}
}
\examples{
\dontrun{
# Load images
GALEX_NUV=readFITS(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))
VST_r=readFITS(system.file("extdata", 'VST_r.fits', package="magicaxis"))
VISTA_K=readFITS(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))

# Warp to common WCS:
GALEX_NUV_VST=magwarp(GALEX_NUV, VST_r$hdr)
VISTA_K_VST=magwarp(VISTA_K, VST_r$hdr)

# Run profoundMultiBand:
multi=profoundMultiBand(inputlist=list(GALEX_NUV_VST, VST_r, VISTA_K_VST),
magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'))

# Check some of the outputs:
plot(multi$pro_detect)

magplot(multi$cat_tot$mag_rt, multi$cat_col$mag_NUVc-multi$cat_col$mag_rc, ylim=c(-2,10))
points(multi$cat_tot$mag_rt, multi$cat_col$mag_rc-multi$cat_col$mag_Kc, col='red')

# Some options on passing segim:

multi2=profoundMultiBand(segim=multi$pro_detect$segim, inputlist=list(GALEX_NUV_VST,
VST_r, VISTA_K_VST), magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'),
iters_det = 0)

multi3=profoundMultiBand(segim=multi$pro_detect$segim_orig, inputlist=list(GALEX_NUV_VST,
VST_r, VISTA_K_VST), magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'),
iters_det = 6)

#multi and multi3 should create identical plots (since we are dilating the original
#segim_orig in the same manner), but multi2 will just be the final dilated segim without
#any dilations, hence the Dilations histogram is all 0. The final fluxes should be the
#same though (bottom-centre and bottom-right).

plot(multi)
plot(multi2)
plot(multi3)
}
}
\keyword{ Detection }% use one of  RShowDoc("KEYWORDS")
\keyword{ Photometry }% __ONLY ONE__ keyword per line
