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
profoundFitMagPSF(xcen = NULL, ycen = NULL, mag = NULL, image = NULL, sigma = NULL,
mask = NULL, psf = NULL, iters = 5, magdiff = 1, modxy = FALSE, sigthresh = 0,
itersub = TRUE, magzero = 0, modelout = TRUE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Numeric vector; required, x centres. Must be the same length as \option{mag}. These values remain fixed during the fitting.
}
  \item{ycen}{
Numeric vector; required, y centres. Must be the same length as \option{mag}. These values remain fixed during the fitting.
}
  \item{mag}{
Numeric vector; required, initial PSF mags. This is what will be fitted.
}
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{sigma}{
Sigma matrix; required, the measurement errors per pixel (expressed in terms of sigma).
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{psf}{
Numeric matrix; required, an empirical point spread function (PSF) image matrix that \code{ProFit} will use to convolve the model.
}
  \item{iters}{
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
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since a big \option{image} can take a long time to run, you might want to monitor progress.  
}
}
\details{
This function uses \code{ProFit} to make a full model image. It then makes a model image for each indiviual PSF component and optimises just this alone. This means a maximum likelihood solution is converged on very efficiently, and is inspired by Expectation Maximisation which is often used for mixture model problems (which this basically is) for the reason of fast convergence. Here we use \code{\link{optim}} with \option{method} Brent to achieve the individual PSF optimisations.
}
\value{
List containing:

\item{psfstats}{A data.frame of photometric properties (see below).}
\item{origmodel}{The original model image before optimisation. This will have the same dimensions as the input \option{image}. Only relevant if \option{modelout}=TRUE, else NA.}
\item{finalmodel}{The final model image after optimisation. This will have the same dimensions as the input \option{image}. Only relevant if \option{modelout}=TRUE, else NA.}
\item{origLL}{The original data-model log-likelihood. Only relevant if \option{modelout}=TRUE, else NA.}
\item{finalLL}{The final data-model log-likelihood. Only relevant if \option{modelout}=TRUE, else NA.}

The \option{psfstats} data.frame has the following columns:

\item{xcen}{Numeric vector; x centres. If \option{modxy}=FALSE these will be the same as the input \option{xcen}, but if \option{modxy}=TRUE they will be adjusted.}
\item{ycen}{Numeric vector; y centres. If \option{modxy}=FALSE these will be the same as the input \option{ycen}, but if \option{modxy}=TRUE they will be adjusted.}
\item{flux}{Numeric vector; the final fitted PSF fluxes. This will be the same length as the input \option{mag} vector.}
\item{flux_err}{Numeric vector; the final fitted PSF flux errors. This will be the same length as the input \option{mag} vector.}
\item{mag}{Numeric vector; the final fitted PSF magnitudes. This will be the same length as the input \option{mag} vector.}
\item{mag_err}{Numeric vector; the final fitted PSF magnitude errors. This will be the same length as the input \option{mag} vector.}
\item{psf}{Numeric vector; the log-likelihood of the individual source fit.}
\item{signif}{Numeric vector; approximate singificance of the detection using the Chi-Square distribution.}
}
\references{
Expectation Maximisation \href{http://en.wikipedia.org/wiki/Expectation-maximization_algorithm}{Wikipedia}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundFluxDeblend}}
}
\examples{
\dontrun{
s250_im=readFITS(system.file("extdata", 'IRdata/s250_im.fits', package="ProFound"))
s250_psf=readFITS(system.file("extdata", 'IRdata/s250_psf.fits', package="ProFound"))
s250_psf_cut=magcutout(s250_psf$imDat, box=c(21,21))$image

pro_s250=profoundProFound(s250_im, pixcut=1, skycut=2, ext=1, redosky=FALSE, iters=1,
tolerance=0, sigma=0)
pro_s250$segstats=pro_s250$segstats[!is.na(pro_s250$segstats$mag),]

newmag=profoundFitMagPSF(xcen=pro_s250$segstats$xcen, ycen=pro_s250$segstats$ycen,
mag=pro_s250$segstats$mag, image=s250_im$imDat-pro_s250$sky, sigma=pro_s250$skyRMS,
psf=s250_psf_cut)

magimage(s250_im$imDat-newmag$origmodel, magmap=FALSE, col=topo.colors(1e3),
zlim=c(-0.02,0.02))
magimage(s250_im$imDat-newmag$finalmodel, magmap=FALSE, col=topo.colors(1e3),
zlim=c(-0.02,0.02))

magplot(pro_s250$segstats$mag, newmag$mag,xlim=c(-1,6),ylim=c(-1,6),xlab='Start Mag',
ylab='End Mag', asp=1, grid=TRUE)
magerr(pro_s250$segstats$mag, newmag$mag, xlo=pro_s250$segstats$mag_err,
ylo=newmag$magerr)
abline(0,1, col='red')

#We can also run slightly adjusting the xcen and ycen as we go:

newmag2=profoundFitMagPSF(xcen=pro_s250$segstats$xcen, ycen=pro_s250$segstats$ycen,
mag=pro_s250$segstats$mag, image=s250_im$imDat-pro_s250$sky, sigma=pro_s250$skyRMS,
psf=s250_psf_cut, modxy=TRUE)

magimage(s250_im$imDat-newmag2$finalmodel, magmap=FALSE, col=topo.colors(1e3),
zlim=c(-0.02,0.02))

magplot(pro_s250$segstats$mag, newmag2$mag,xlim=c(-1,6),ylim=c(-1,6),xlab='Start Mag',
ylab='End Mag', asp=1, grid=TRUE)
magerr(pro_s250$segstats$mag, newmag$mag, xlo=pro_s250$segstats$mag_err,
ylo=newmag2$magerr)
abline(0,1, col='red')

# The two profoundFitMagPSF approaches agree well within the error:

magplot(newmag$mag, newmag2$mag,xlim=c(-1,6),ylim=c(-1,6),xlab='Start Mag',
ylab='End Mag', asp=1, grid=TRUE)
magerr(newmag$mag, newmag2$mag, xlo=newmag$magerr, ylo=newmag2$magerr)
abline(0,1, col='red')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ deblend }% use one of  RShowDoc("KEYWORDS")
