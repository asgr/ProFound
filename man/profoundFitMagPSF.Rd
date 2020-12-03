\name{profoundFitMagPSF}
\alias{profoundFitMagPSF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit PSF Magnitudes to an Image
}
\description{
Fits PSF mags to an image with known source positions. This is the best deblend option in the regime where the sources are not well resolved, but the source positions are reasonably well defined.
}
\usage{
profoundFitMagPSF(xcen = NULL, ycen = NULL, RAcen = NULL, Deccen = NULL, mag = NULL,
image = NULL, im_sigma = NULL, mask = NULL, psf = NULL, fit_iters = 5, magdiff = 1,
modxy = FALSE, sigthresh = 0, itersub = TRUE, magzero = 0, modelout = TRUE,
fluxtype = 'Raw', psf_redosky = FALSE, fluxext = FALSE, header = NULL, doProFound = FALSE,
findextra = FALSE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Numeric vector; x centres. If provided, must be the same length as \option{mag}, and be paired with \option{ycen}.
}
  \item{ycen}{
Numeric vector; y centres. If provided, must be the same length as \option{mag}, and be paired with \option{xcen}.
}
  \item{RAcen}{
Numeric vector; right ascension centres in degrees.  If provided, must be the same length as \option{mag}, and be paired with \option{Deccen}.
}
  \item{Deccen}{
Numeric vector; declination centres in degrees.  If provided, must be the same length as \option{mag}, and be paired with \option{RAcen}.
}
  \item{mag}{
Numeric vector; required, initial PSF mags. This is what will be fitted.
}
  \item{image}{
Numeric matrix; required, the image we want to analyse. If \option{image} is a list as created by \code{readFITS}, \code{read.fits} of \code{\link{magcutoutWCS}} then the image part of these lists is passed to \option{image} and the correct header part is passed to \option{header}. Note, image NAs are treated as masked pixels.
}
  \item{im_sigma}{
Numeric matrix; required, the measurement errors per pixel (expressed in terms of sigma).
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{psf}{
Numeric matrix; required, an empirical point spread function (PSF) image matrix that \code{ProFit} will use to convolve the model.
}
  \item{fit_iters}{
Integer scalar; how many iterations should be run? Usually converges quite quickly, so rarely needs to be much larger than 5 (default).
}
  \item{magdiff}{
Numeric scalar; required, what is the allowed magnitude adjustment per iteration. Smaller values means faster fitting, but the correct solution obviously needs to lie with \option{iters}x\option{magdiff} for all PSF magnitudes.
}
  \item{modxy}{
Logical; should \option{xcen} and \option{ycen} positions be adjusted during fitting? This adds computation time (roughly an extra 50\%), but generally produces slightly better fits. It is limited so positions can only move 0.5 pixels per iteration, so the maximum move possible is \option{iters} x 0.5 pixels.
}
  \item{sigthresh}{
Numeric scalar; the cut level to apply before re-estimating the \option{xcen} and \option{ycen} positions. Higher means only brighter pixels are considered so is more robust, but it means fainter sources might not be adjusted at all. Only relevant if \option{modxy}=TRUE.
}
  \item{itersub}{
Logical; should each marginalised source profile be subtracted as we loop around the sources within a given iteration? The reason to perhaps set this to FALSE is that the source order will then affect the results (the earlier source will collect more flux where they overlap), but the overall solution will always be better with this set to TRUE (default).
}
  \item{magzero}{
Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega).  
}
  \item{modelout}{
Logical; should the full model image be output?
}
  \item{fluxtype}{
Character scaler; specifies whether fluxes will be output in Jansky / MicroJansky ('Jansky' / 'microjansky'), or in raw/unscaled image ADUs ('Raw' / 'ADU' / 'ADUs', the default). You can only use 'Jansky' if the specified \option{magzero} gets the data into the AB system, else the fluxes will not be Jansky.
}
  \item{psf_redosky}{
Logical; should the sky calulcated by \code{\link{profoundProFound}} be subtracted from the target image? Since FIR frames are often confused, this should be turned on with care since possibly none of the pixels are a good sky reference. Put another way, you really need to trust your PSF wings.
}
  \item{fluxext}{
Logical; should the full model be used to compute an image deblending function better suited for extended sources? This means fluxes are guaranteed to add up to those available in the \option{image}. The fit PSF can differ, especially if the \option{psf} supplied is not in detail identical to the true \option{image} PSF. In practice we find that using \option{fluxext}=FALSE works better when the sources are not well resolved (assuming the true \option{psf} ~ true PSF), but when sources are marginally extended then \option{fluxext}=TRUE is probably the safer mode to use (hence the name).
}
  \item{header}{
Full FITS header in table or vector format. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the \code{astro} package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. If a header is provided then key words will be taken from here as a priority. Missing header keywords are printed out and other header option arguments are used in these cases.
}
  \item{doProFound}{
Logical; if TRUE then \code{\link{profoundProFound}} will be run with the arguments provided by \dots. This will be used to determine \option{xcen}, \option{ycen} and \option{mag} blindly.
}
  \item{findextra}{
Logical; if TRUE then after an initial run of \code{profoundFitMagPSF} \code{\link{profoundProFound}} is used to find additional sources that cannot be well modelled. The combined set if the initial sources and the new extra sources are then passed back into \code{profoundFitMagPSF} for a final run. The main source list is output to \option{psfstats}, and the extra sources are output to \option{psfstats_extra}. \option{psfstats_extra} is otherwise NULL.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since a big \option{image} can take a long time to run, you might want to monitor progress.  
}
  \item{\dots}{
Arguments to be passed on to \code{\link{profoundProFound}}. Only relevant if \option{doProFound}=TRUE or \option{findextra}=TRUE.
}
}
\details{
This function uses \code{ProFit} to make a full model image. It then makes a model image for each indiviual PSF component and optimises just this alone. This means a maximum likelihood solution is converged on very efficiently, and is inspired by Expectation Maximisation which is often used for mixture model problems (which this basically is) for the reason of fast convergence. Here we use \code{\link{optim}} with \option{method} Brent to achieve the individual PSF optimisations.
}
\value{
Object of class "fitmagpsf", a list containing:

\item{psfstats}{Data.frame; main source photometric properties (see below).}
\item{origmodel}{Numeric matric; the original model image before optimisation. This will have the same dimensions as the input \option{image}. Only relevant if \option{modelout}=TRUE, else NA.}
\item{finalmodel}{Numeric matric; the final model image after optimisation. This will have the same dimensions as the input \option{image}. Only relevant if \option{modelout}=TRUE, else NA.}
\item{origLL}{Numeric scalar; the original data-model log-likelihood. Only relevant if \option{modelout}=TRUE, else NA.}
\item{finalLL}{Numeric scalar; the final data-model log-likelihood. Only relevant if \option{modelout}=TRUE, else NA.}
\item{image}{Numeric matric; the image directly used for fitting. This might have an additional sky component removed (if \option{doProFound}=TRUE) and NAs for masked regions.}
\item{header}{The header provided, if missing this is NULL.}
\item{psfstats_extra}{Data.frame; xtra photometric properties if \option{findextra}=TRUE, otherwise NULL (see below).}
\item{profound}{List; an object of class profound as output by the \code{\link{profoundProFound}} function. This contains the output for the ProFound part of the analysis if \option{doProFound}=TRUE or \option{findextra}=TRUE.}
\item{mask}{The input or computed \option{mask} matrix, else NULL.}
\item{call}{The original function call.}
\item{date}{The date, more specifically the output of \code{\link{date}}.}
\item{time}{The elapsed run time in seconds.}
\item{ProFound.version}{The version of \code{ProFound} run, more specifically the output of \code{\link{packageVersion}}('ProFound').}
\item{R.version}{The version of \code{R} run, more specifically the output of \code{\link{R.version}}.}
  

\option{psfstats} \option{psfstats_extra} data.frames have the following columns:

\item{xcen}{Numeric vector; x centres. If \option{modxy}=FALSE these will be the same as the input \option{xcen}, but if \option{modxy}=TRUE they will be adjusted.}
\item{ycen}{Numeric vector; y centres. If \option{modxy}=FALSE these will be the same as the input \option{ycen}, but if \option{modxy}=TRUE they will be adjusted.}
\item{flux}{Numeric vector; the final fitted PSF fluxes (either Jansky or raw, depending on \option{fluxtype}). This will be the same length as the input \option{mag} vector.}
\item{flux_err}{Numeric vector; the final fitted PSF flux errors (either Jansky or raw, depending on \option{fluxtype}). This will be the same length as the input \option{mag} vector.}
\item{mag}{Numeric vector; the final fitted PSF magnitudes. This will be the same length as the input \option{mag} vector.}
\item{mag_err}{Numeric vector; the final fitted PSF magnitude errors. This will be the same length as the input \option{mag} vector.}
\item{psf}{Numeric vector; the log-likelihood of the individual source fit.}
\item{signif}{Numeric vector; approximate singificance of the detection using the Chi-Square distribution.}
}
\references{
Expectation Maximisation \href{https://en.wikipedia.org/wiki/Expectation-maximization_algorithm}{Wikipedia}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundFluxDeblend}}, \code{\link{plot.fitmagpsf}}
}
\examples{
\dontrun{
s250_im=readFITS(system.file("extdata", 'IRdata/s250_im.fits', package="ProFound"))
s250_psf=readFITS(system.file("extdata",'IRdata/s250_psf.fits', package="ProFound"))$imDat
magzero_s250=11.68

pro_s250=profoundProFound(s250_im, pixcut=1, skycut=2, ext=1, redosky=FALSE, iters=1,
tolerance=0, sigma=0, magzero=magzero_s250)
pro_s250$segstats=pro_s250$segstats[!is.na(pro_s250$segstats$mag),]

newmag=profoundFitMagPSF(RAcen=pro_s250$segstats$RAcen, Deccen=pro_s250$segstats$Deccen,
image=s250_im, psf=s250_psf, doProFound=TRUE, findextra=TRUE, verbose=TRUE, redosky=FALSE,
magzero=magzero_s250)

magimage(newmag$image, qdiff=TRUE)
magimage(newmag$image - newmag$origmodel, qdiff=TRUE)
magimage(newmag$image - newmag$finalmodel, qdiff=TRUE)

magplot(pro_s250$segstats$mag, newmag$psfstats$mag, xlim=c(11,18), ylim=c(11,18),
xlab='ProFound Mag', ylab='FitPSF Mag', asp=1, grid=TRUE)
magerr(pro_s250$segstats$mag, newmag$psfstats$mag, xlo=pro_s250$segstats$mag_err,
ylo=newmag$psfstats$mag_err)
abline(0,1, col='red')

#We can also run slightly adjusting the xcen and ycen as we go:

newmag2=profoundFitMagPSF(RAcen=pro_s250$segstats$RAcen, Deccen=pro_s250$segstats$Deccen,
image=s250_im, psf=s250_psf, doProFound=TRUE, findextra=TRUE, verbose=TRUE, redosky=FALSE,
modxy=TRUE, magzero=magzero_s250)

magimage(newmag2$image - newmag2$finalmodel, qdiff=TRUE)

magplot(pro_s250$segstats$mag, newmag2$psfstats$mag,xlim=c(11,18), ylim=c(11,18),
xlab='ProFound Mag', ylab='FitPSF Mag', asp=1, grid=TRUE)
magerr(pro_s250$segstats$mag, newmag2$psfstats$mag, xlo=pro_s250$segstats$mag_err,
ylo=newmag2$psfstats$mag_err)
abline(0,1, col='red')

# The two profoundFitMagPSF approaches agree well within the error:

magplot(newmag$psfstats$mag, newmag2$psfstats$mag,xlim=c(11,18), ylim=c(11,18),
xlab='FitPSF Mag', ylab='FitPSF (mod xy) Mag', asp=1, grid=TRUE)
magerr(newmag$psfstats$mag, newmag2$psfstats$mag, xlo=newmag$psfstats$mag_err,
ylo=newmag2$psfstats$mag_err)
abline(0,1, col='red')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ deblend }% use one of  RShowDoc("KEYWORDS")
