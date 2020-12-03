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
profoundMultiBand(inputlist = NULL, dir = "", segim = NULL, mask = NULL,
detectbands = "r", multibands = c("u", "g", "r", "i", "z"), iters_det = 6, iters_tot = 0,
sizes_tot = 5, magzero = 0, gain = NULL, box = 100, grid=box, boxadd=box/2, app_diam = 1,
bandappend = multibands, totappend = "t", colappend = "c", grpappend = 'g', dotot = TRUE,
docol = TRUE, dogrp = TRUE, deblend = FALSE, groupstats = FALSE,
groupby_det = 'segim_orig', groupby_mul = 'segim_orig', keepsegims = FALSE,
masking = 'and', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputlist}{
A list of already loaded images. Typically of the type loaded in from FITS files by the \code{astro} package's \code{read.fits} function, or the \code{FITSio} package's \code{readFITS} function. If using the \option{inputlist} parameter the length of the list must be the same length as \option{multibands} (and the related parameters).
}
  \item{dir}{
If \option{inputlist} is left as NULL then \code{profoundMultiBand} will instead try to load in FITS images from the directory specified by \option{dir}. The images in the directory must have names like \option{multibands}[i].fits etc (so with the defaults names like u.fits and g.fits would be okay). Since \option{multibands} effectively specifies the file names much more complicated naming can be used and passed in, but it is also used by default for naming the catalogue column outputs, so shorter names/references are likely to be preferable there (i.e. mag_ut is simpler than mag_KiDS_VST_ut etc). This can be over-ridden by using \option{bandappend}.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This matrix *must* be the same dimensions as the detection image/s if supplied. If this option is used then \code{profoundMultiBand} will not compute its initial segmentation map using \code{\link{profoundMakeSegim}}, which is then dilated. Instead it will use the one passed through \option{segim} and dilate this according to the \option{iters_det} argument (so set this to 0 if you want the \option{segim} to be used as is).
}
  \item{mask}{
Boolean matrix or integer scalar; optional, parts of the image to mask out (i.e. ignore). If a matrix is provided, this matrix *must* be the same dimensions as \option{image} where 1 means mask out and 0 means use for analysis. if a scalar is provided it indicates the exact \option{image} values that should be treated as masked (e.g. by setting masked pixels to 0 or -999). The latter achieves the same effect as setting masked \option{image} pixels to NA, but allows for the fact not all programs can produce \code{R} legal NA values. Not that for detection masking is an OR operator, so if a pixel is masked on any input \option{image} it will be masked in the detection \option{image}. This is the more conservative way to treat the issue of multiband masking.
}
  \item{detectbands}{
Character vector; the names of the detection bands that will be stacked using \code{\link{profoundMakeStack}} and then analysed with the provided settings with \code{\link{profoundProFound}} to make a reference segmentation map for further multi band photometry. These bands must be present in \option{multibands}. Can be a scalar (i.e. a single band is used). If set to 'get' then it will use all legal FITS files in the target directory. If set to 'all' then it will use all \option{multibands} inputs.
}
  \item{multibands}{
Character vector; the names of the target multi band photometry images. If set to 'get' then it will use all legal FITS files in the target directory. If using the \option{inputlist} parameter the length of the list must be the same length as \option{multibands}. \option{magzero} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. If specified, \option{gain} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. If specified, \option{catappend} must either the same length and order as \option{multibands}.
}
  \item{iters_det}{
Integer scalar; the maximum number of curve of growth dilations that should be made to the detection image. This needs to be large enough to capture all the flux for sources of interest, but increasing this will increase the computation time for \code{profoundProFound}. If this is set to 0 then the undilated \option{segim} image, whether provided or computed internally via \code{\link{profoundMakeSegim}}, will be used instead.
}
  \item{iters_tot}{
Integer vector; the maximum number of curve of additional growth dilations that should be made above the dilated detection segmentation map for multi band total colour photometry. This is only relevant if \option{dotot}=TRUE. This should not be set too high (and might even be 0, the default) since the detection image should generally be fairly deep. \option{iters_tot} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}.
}
  \item{sizes_tot}{
Integer vector;  the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number. See \code{makeBrush}. Passed to \code{\link{profoundMakeSegimDilate}}. \option{sizes_tot} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}.
}
  \item{magzero}{
Numeric vector; the magnitude zero point of the images being used. \option{magzero} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{gain}{
Numeric vector; the gain of the images being used. \option{gain} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{box}{
Numeric vector; the sky estimate box size of the images being used. \option{box} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{grid}{
Integer vector; the resolution of the background grid to estimate the sky with. By default this is set to be the same as the \option{box}. \option{grid} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{boxadd}{
Integer vector; the dimensions to add to the \option{box} to capture more pixels if \option{skypixmin} has not been achieved. By default this is set to be the same as the \option{box}/2. \option{boxadd} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
}
  \item{app_diam}{
Numeric vector; the diameter in arc seconds to use for pseudo aperture photometry. This will use the appropriate pixel scale to convert the aperture into image units. The psuedo aperture photometry is output to columns \option{flux_app} and \option{mag_app} in \option{cat_tot}. \option{app_diam} must either be length 1 (in which this value is used for all bands), or the same length and order as \option{multibands}. See also \code{\link{profoundProFound}}.
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
  \item{deblend}{
Logical; should total segment flux be deblended using \code{\link{profoundFluxDeblend}} and these columns appended to the end of the output segstats? This only applies to the \option{cat_tot} output.
}
  \item{groupstats}{
Logical; if TRUE then the IDs of grouped segments is calculated for the detection image via \code{\link{profoundSegimGroup}} and output to the returned object \option{group}. By default this option is linked to \option{boundstats}, i.e. it is assumed if you want boundary statistics then you probably also want grouped object IDs returned.
}
  \item{groupby_det}{
Character scalar; How should the grouped segmentation map be formed that will be used to produce the \option{cat_grp} output? Options are either via 'segim' or 'segim_orig'. 'segim' will create more groups, 'segim_orig' will have less.
}
  \item{groupby_mul}{
Character scalar; How should the grouped segmentation map be formed that will be used to produce the multi-band deblending for \option{cat_tot} and \option{cat_col} outputs? Options are either via 'segim' or 'segim_orig'. 'segim' will create more groups, 'segim_orig' will have less.
}
  \item{keepsegims}{
Logical; if TRUE then \code{profoundMultiBand} will return a list of per band segmentation maps in \option{segimlist}. This is useful since they can be slightly different per band.
}
  \item{masking}{
Character scalar; what to do with masked pixels (NAs in the \option{image}). If 'or' a pixel is masked if *any* of the images being stacked have a masked pixel (NA in the \option{image_list}) at that location, if 'and' then a pixel is masked if *all* of the images being stacked have a masked pixel at that location.
}
  \item{\dots}{
Further arguments to be passed to all instances of \code{\link{profoundProFound}}. E.g. if the sky \option{SBdilate} is set to a value (default is NULL) this will be propogated to all of the multi band photometry runs of \code{\link{profoundProFound}}.
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
\item{segimlist}{If \option{keepsegims}=TRUE then this object contains a list of per band total segmentation maps. This is useful since they can be slightly different per band.}
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

# Run profoundMultiBand on defaults:
multi=profoundMultiBand(inputlist=list(GALEX_NUV_VST, VST_r, VISTA_K_VST),
magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'))

# Notice the blue halo around the central sources:
plot(multi$pro_detect)

# Run profoundMultiBand with boxiters=2 (to avoid over-subtracting the sky):
multi=profoundMultiBand(inputlist=list(GALEX_NUV_VST, VST_r, VISTA_K_VST),
magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'), boxiters = 2)

# Looks better now:
plot(multi$pro_detect)

magplot(multi$cat_tot$mag_rt, multi$cat_col$mag_NUVc-multi$cat_col$mag_rc, ylim=c(-2,10))
points(multi$cat_tot$mag_rt, multi$cat_col$mag_rc-multi$cat_col$mag_Kc, col='red')

# Some options on passing segim:

multi2=profoundMultiBand(segim=multi$pro_detect$segim, inputlist=list(GALEX_NUV_VST,
VST_r, VISTA_K_VST), magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'),
iters_det = 0, boxiters=2)

multi3=profoundMultiBand(segim=multi$pro_detect$segim_orig, inputlist=list(GALEX_NUV_VST,
VST_r, VISTA_K_VST), magzero=c(20.08,0,30), detectbands='r', multibands=c('NUV','r','K'),
iters_det = 6, boxiters=2)

# multi and multi3 should create identical plots (since we are dilating the original
# segim_orig in the same manner), but multi2 will just be the final dilated segim without
# any dilations, hence the top-right is all green (segim=segim_orig). The final fluxes
# should be the same though for all 3 runs (left-middle, bottom-centre and bottom-right).

plot(multi$pro_detect)
plot(multi2$pro_detect)
plot(multi3$pro_detect)
}
}
\concept{ Detection }% use one of  RShowDoc("KEYWORDS")
\concept{ Photometry }% __ONLY ONE__ keyword per line
