\name{profoundAperPhot}
\alias{profoundAperPhot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aperture Fluxes
}
\description{
Allows you to extract one or many aperture fluxes from a target image.
}
\usage{
profoundAperPhot(image = NULL, segim = NULL, app_diam = 1, tar = NULL,
  pixscale = 1, magzero = 0, correction = TRUE, verbose = FALSE)
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
  \item{verbose}{
Logical; should verbose output be displayed to the user?
}
}
\details{
This reproduces the output of the \code{\link{profoundProFound}} \option{app_diam} aperture magnitude, but the advantage is being able to specify multiple \option{app_diam} and also target only certain segID and specify the precise xcen and ycen you wish to use.

Note by default \option{flux_app_X} and \option{mag_app_X} are already corrected by the fractional coverage \option{frac_app_X}. This is done by assessing how many pixels over or under the desired area the pixelated aperture is, and then using the mean of the outer aperture pixel flux to adjust \option{flux_out} (either increasing the raw flux if \option{frac_app_X} < 1, or decreasing the raw flux if \option{frac_app_X} > 1). For this reason users should not attempt their own additional correction when \option{correction} = TRUE.
}
\value{
Data.frame containing:

\item{segID}{Segmentation ID}
\item{flux_app_X}{The raw flux (flux_app_1 etc)} 
\item{mag_app_X}{The magnitude (mag_app_1 etc)}
\item{N_app_X}{The number of pixels falling in the aperture (N_app_1 etc)}
\item{frac_app_X}{The area of \option{N_app_X} divided by target aperture area (N_app_1 etc)}
\item{flux_min}{The mean pixel flux of the outer most (usually faintest) pixels in the aperture}

Note the '_X' refers to the aperture by position in the input \option{app_diam} vector.
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
  app_diam=1:3, magzero=30)

#These should look the same
magplot(profound$segstats$mag_app, profound$segstats$mag_app -new_app$mag_app_2,
  ylim=c(-0.1,0.1), xlab='Pro mag_app', ylab='Pro - Aper (mag_app)')

tar = data.frame(segID=1:5)
new_app2 = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
  app_diam=1:3, magzero=30, tar=tar)
new_app2
  
#we can also target specific regions (but these then need to overlap with the segID):

tar2 = profound$segstats[1:8, c("segID", "xcen", "ycen")]
new_app3 = profoundAperPhot(image=image - profound$sky, segim=profound$segim,
                            app_diam=1:3, magzero=30, tar=tar2)
new_app3
}
\concept{ aperture }
