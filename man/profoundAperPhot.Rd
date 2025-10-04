\name{profoundAper}
\alias{profoundAperPhot}
\alias{profoundAperRan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aperture Fluxes
}
\description{
Allows you to extract one or many aperture fluxes from a target image.
}
\usage{
profoundAperPhot(image = NULL, segim = NULL, app_diam = 1, keyvalues = NULL,
  tar = NULL, pixscale = 1, magzero = 0, correction = TRUE, centype = 'mean',
  fluxtype = 'Raw', verbose = FALSE)
  
profoundAperRan(image = NULL, segim = NULL, app_diam = 1, Nran = 100,
  keyvalues = NULL, pixscale = 1, magzero = 0, correction = TRUE,
  centype = 'mean', fluxtype = 'Raw', verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. If \option{image} is a list as created by \code{Rfits_read_image}  then the image part of these lists is passed to \option{image} and the correct keyvalues part is passed to \option{keyvalues}. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This is required.
}
  \item{app_diam}{
Numeric vector; the diameter in arc seconds to use for circular aperture photometry (hence 'app'). This will use the appropriate pixel scale to convert the aperture into image units. The circular aperture photometry is output to columns \option{flux_app} and \option{mag_app} appended by the vector number.
}
  \item{Nran}{
Integer scalar; the number of random apertures to place in the sky (where \option{segim} is 0).  
}
  \item{keyvalues}{
List; header values to be used for the WCS. This can be useful if you have split up your matrix image data (passing that to \option{image}) but still have standalone \option{keyvalues} available to use. Note if \option{keyvalues} are found attached to the \option{image}
}
  \item{tar}{
Data.frame; allows you to target only specific \option{segID} and optionally provide the target \option{xcen} and \option{ycen} to use (otherwise these are determined automatically based on the brightest pixel in each segment). Note if providing positions both \option{xcen} and \option{ycen} must be provided for all segIDs (nothing should be set to NA). See Examples.
}
  \item{pixscale}{
Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). Note if this is not set then an attempt is made to determine the \option{pixelscale} from the detected WCS information.
}
  \item{magzero}{
Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega).
}
  \item{correction}{
Logical; should discreteness and fractional coverage corrections be applied to the aperture flux (changing the output \option{flux_app_X} and \option{mag_app_X} columns). Unless you really know what you are doing, this should be left on TRUE (and in this case do not apply additional corrections). See Details.
}
  \item{centype}{
Character scalar; type of automatic centre to use per segment. Either 'max' (the location of the peak pixel flux in the segment), or 'mean' (the weighted flux centre). 'max' will always be exactly in the centre of the peak flux pixel, whereas 'mean' can be any fractional pixel position. This can impact the aperture photometry in subtle ways. 
}
  \item{fluxtype}{
Character scaler; specifies whether fluxes will be output in Jansky / MicroJansky ('Jansky' / 'microjansky'), or in raw/unscaled image ADUs ('Raw' / 'ADU' / 'ADUs', the default). You can only use Jansky / MicroJansky if the specified \option{magzero} gets the data into the AB system, else the fluxes will not be Jansky.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user?
}
}
\details{
\code{profoundAperPhot} reproduces the output of the \code{\link{profoundProFound}} \option{app_diam} aperture magnitude, but the advantage is being able to specify multiple \option{app_diam} and also target only certain segID and specify the precise xcen and ycen you wish to use.

Note by default \option{flux_app_X} and \option{mag_app_X} are already corrected by the fractional coverage \option{frac_app_X}. This is done by assessing how many pixels over or under the desired area the pixelated aperture is, and then using the mean of the outer aperture pixel flux to adjust \option{flux_out} (either increasing the raw flux if \option{frac_app_X} < 1, or decreasing the raw flux if \option{frac_app_X} > 1). For this reason users should not attempt their own additional correction when \option{correction} = TRUE.

\code{profoundAperRan} places random apertures on the sky where \option{segim} is 0. Any non-zero pixels in \option{segim} are forced to be zero in the created \option{segim_ran}. It then returns the output of \code{profoundAperPhot} at these aperture positions, passing all the additional matching argument on. This can be used to asses pragmatic aperture errors based on the \option{image}, accounting for sky variability, skyRMS, undetected background sources, and covariance in pixel flux. The computed \option{errors} uses only full apertures (so no pixels lost to neighbour coverage). Depending on your use case, you might want to only consider using apertures with partial coverage.
}
\value{
\code{profoundAperPhot} returns data.frame containing:

\item{segID}{Segmentation ID}
\item{xcen}{The aperture x position}
\item{ycen}{The aperture y position}
\item{flux_app_X}{The raw flux (flux_app_1 etc)} 
\item{mag_app_X}{The magnitude (mag_app_1 etc)}
\item{SB_app_X}{The surface brightness in mag/asec^2 units (SB_app_1 etc)}
\item{N_app_X}{The number of pixels falling in the aperture (N_app_1 etc)}
\item{frac_app_X}{The area of \option{N_app_X} divided by target aperture area (N_app_1 etc)}
\item{flux_min}{The mean pixel flux of the outer most (usually faintest) pixels in the aperture}

Note the '_X' refers to the aperture by position in the input \option{app_diam} vector.

\code{profoundAperRan} returns a list containing:

\item{AperPhot}{The full output of \code{profoundAperPhot} run on the random apertures (in case you want to compute different error statistics compared to \option{errors}.}
\item{errors}{The 16-84 quantile based error in aperture fluxes, in units of \option{fluxtype}. The same length and order as \option{app_diam}.}
\item{segim_ran}{The random apertures used, placed on the segmentation map.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}
}
\examples{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))

profound=profoundProFound(image, magzero=30, app_diam=2)

new_app = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
  app_diam=1:3, magzero=30, verbose=TRUE)

#These should look the same
magplot(profound$segstats$mag_app, profound$segstats$mag_app -new_app$mag_app_2,
  ylim=c(-0.1,0.1), xlab='Pro mag_app', ylab='Pro - Aper (mag_app)')

tar2 = data.frame(segID=1:5)
new_app2 = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
  app_diam=1:3, magzero=30, tar=tar2, verbose=TRUE)
new_app2
  
#we can also target specific regions (but these then need to overlap with the segID).
#here we use xcen and ycen, which is different the default centype (max):

tar3 = profound$segstats[1:8, c("segID", "xcen", "ycen")]
new_app3 = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
                            app_diam=1:3, magzero=30, tar=tar3, verbose=TRUE)
new_app3

#here we use xmax and ymax, which is the same as the default centype (max):

tar4 = profound$segstats[1:8, c("segID", "xmax", "ymax")]
new_app4 = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
                            app_diam=1:3, magzero=30, tar=tar4, verbose=TRUE)
new_app4

#Different positions do make some difference:
magplot(new_app3$mag_app_2, new_app3$mag_app_2 - new_app4$mag_app_2,
  ylim=c(-0.1,0.1), xlab='Pro mag_app', ylab='Cen - Max (mag_app)')
  
# Aperture errors:
new_error = profoundAperRan(image=image - profound$sky, segim=profound$segim,
  app_diam=1:3, magzero=30, verbose=TRUE)
new_error$errors

magplot(profound$segstats$mag, new_app$flux_app_1, log='y', ylim=c(10,2e6),
  xlab='mag', ylab='flux_app')
magerr(profound$segstats$mag, new_app$flux_app_1, ylo=new_error$errors[1])

#mag errors will be given by:
maghist((2.5/log(10))*new_error$errors[1]/new_app$flux_app_1,
  xlab='Error [mag]')
}
\concept{ aperture }
