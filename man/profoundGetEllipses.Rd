\name{profoundGetEllipses}
\alias{profoundGetEllipses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate multiple annulus properties of iso-photal ellipses
}
\description{
Returns multiple ellipse properties for an image, assumed to be monotonically decreasing in flux from a bright centre (i.e. a classic galaxy).
}
\usage{
profoundGetEllipses(image = NULL, segim = NULL, segID = 1L, levels = 10, magzero = 0,
  pixscale = 1, fixcen = TRUE, dobox = FALSE, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{segim}{
Integer matrix; optional, the segmentation map of the image. This matrix *must* be the same dimensions as \option{image}.
}
  \item{segID}{
Integer scalar; optional, the desired \option{segim} segment to extract from the \option{image}.
}
  \item{levels}{
Integer scalar or vector. If a scalar this is the number of ellipse levels to extract from the \option{image}. If a vector this specfies the extremes of all fractional levels, i.e. it should generally start at 0 and end at 1 to capture all isophotal levels.
}
  \item{magzero}{
Numeric scalar; the magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega). If provided along with \option{pixscale} then the surface brightness output will represent mag/asec^2.
}
  \item{pixscale}{
Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1 (default), then the output \option{radhi}, \option{radlo} and \option{radav} is in terms of pixels, otherwise they are in arcseconds. If provided along with \option{magzero} then the surface brightness output will represent mag/asec^2.
}
  \item{fixcen}{
Logical; should the ellipse centres be fixed to a common flux weighted centre?
}
  \item{dobox}{
Logical; should boxiness be computed? If FALSE then boxiness is fixed to be 0. If TRUE then boxiness is computed (and other parameters are refined) using a maximum likelihood method. This is more expensive to compute, so the default is FALSE.
}
  \item{plot}{
Logical; should a diagnostic plot be generated? This plot is generated by the \code{\link{profoundGetEllipsesPlot}} function.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{profoundGetEllipsesPlot}}. Only relevant if \option{plot}=TRUE.
}
}
\details{
This higher level function provides an easy way to extract iso-photal ellipses from an image of a galaxy. How it works somewhat replicates IRAF's ellipse, but it is really present to offer useful initial guesses for bulge and disk geometric properties. It certainly does not guarantee to return the same solution as IRAF (in fact I am not exactly aware of how IRAF computes its ellipses).

Internally it works by rank ordering the pixels of the galaxy and dividing these into equi-spaced quantiles of flux (so each annulus will approximately sum to the same amount of flux). This means that the error for each ellipse will be approximately constant. For each annulus it then runs \code{\link{profoundGetEllipse}} to compute the ellipse properties of what is assumed to be a fairly narrow annulus of pixels. The implicit assumption is that the galaxy flux more-or-less monotonically decreases from the centre, and dividing pixels like this will assure the extraction of common iso-photal ellipses. This assumption works well within the inner 90\% of a galaxy's flux, but isophotes can be quite noisy once the galaxy flux gets close to the sky RMS level. This said, the ellipse returned will on average make sense, and ellipses tend to overlap only in very extreme cases (where the geometry is highly non-elliptical or there are close contaminants).
}
\value{
A list containing:

\item{ellipses}{A data.frame of ellipse properties ordered by radius (see below).}
\item{segellipses}{Integer matrix; the ellipse-wise segmentation map matched pixel by pixel to \option{image}. This allows you to see which specific pixels used to compute each ellipse annulus in \option{ellipses}, where the number in the segmentation map refers to \option{segellipseID}.}

\option{ellipses} is a data.frame of ellipse properties ordered by radius. It has the following columns

\item{segellipseID}{The ellipse segment ID that refers to the segmentation map \option{segellipses}.}
\item{fluxfrac}{The approximate fraction of galaxy flux contained within this ellipse.}
\item{xcen}{The flux weighted x centre of the ellipse.}
\item{ycen}{The flux weighted y centre of the ellipse.}
\item{radhi}{The major axis extent of the ellipse (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)..}
\item{radlo}{The minor axis extent of the ellipse (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)..}
\item{radav}{The average radius of the ellipse (units of \option{pixscale}, so if \option{pixscale} represents the standard asec/pix this will be asec)..}
\item{axrat}{The axial ratio of the ellipse as given by \option{radlo}/\option{radhi}.}
\item{ang}{The angle of the ellipse in the usual ProFit sense, see \code{profitMakeModel}.}
\item{box}{The boxiness of the ellipse in the usual ProFit sense, see \code{profitMakeModel}.}
\item{xsd}{The flux weighted standard deviation in x (always in units of pix).}
\item{ysd}{The flux weighted standard deviation in y (always in units of pix).}
\item{covxy}{The flux weighted covariance in xy (always in units of pix).}
\item{corxy}{The flux weighted correlation in xy (always in units of pix).}
\item{flux}{The flux contained in the segmented pixels associated with this ellipse.}
\item{N}{The number of segmented pixels associated with this ellipse.}
\item{SB}{The mean surface brightness of the pixels associated with this ellipse (if \option{pixscale} has been set correctly then this column will represent mag/asec^2, otherwise it will be mag/pix^2).}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundGetEllipsesPlot}}, \code{\link{profoundGetEllipse}}, \code{\link{profoundDrawEllipse}}
}
\examples{
\dontrun{
# We need the ProFit library to show the profile: library(ProFit)
image = readFITS(system.file("extdata", 'KiDS/G278109fitim.fits',
package="ProFit"))$imDat
segim = readFITS(system.file("extdata", 'KiDS/G278109segim.fits',
package="ProFit"))$imDat
ellipses_nobox = profoundGetEllipses(image=image, segim=segim, levels=20, dobox=FALSE,
pixscale=0.2)
ellipses_box = profoundGetEllipses(image=image, segim=segim, levels=20, dobox=TRUE,
pixscale=0.2)

magplot(ellipses_box$ellipses$radhi[4:19], ellipses_nobox$ellipses$SB[4:19],
ylim=c(25,17), grid=TRUE, type='l')
points(ellipses_box$ellipses$radhi[4:19],ellipses_box$ellipses$SB[4:19])
#A rough bulge+disk surface brightness profile (mean axrat~0.6):
rlocs=seq(1,30,by=0.1)
bulge=profitRadialSersic(rlocs, mag=18.2, re=1.7, nser=3)
disk=profitRadialSersic(rlocs, mag=18, re=13, nser=0.7)
lines(rlocs, profoundFlux2SB(bulge, pixscale=0.2), col='red')
lines(rlocs, profoundFlux2SB(disk, pixscale=0.2), col='blue')
lines(rlocs, profoundFlux2SB(bulge+disk, pixscale=0.2), col='green')
#To get correct magnitudes you would need to modify the components by the axrat
#and pixel scale.

#We can do a better 1D fit with ease:
#Since the ellipses are divided by equi-flux we can minimise sum-square of the SB diff:
sumsq1D=function(par=c(17.6, log10(1.7), log10(3), 17.4, log10(13), log10(0.7)),
rad, SB, pixscale=1){
  bulge=profitRadialSersic(rad, mag=par[1], re=10^par[2], nser=10^par[3])
  disk=profitRadialSersic(rad, mag=par[4], re=10^par[5], nser=10^par[6])
  total=profoundFlux2SB(bulge+disk, pixscale=pixscale)
  return=sum((total-SB)^2)
}

lower=c(10,0,-0.5,10,0,-0.5)
upper=c(30,2,1,30,2,1)

fit1D=optim(sumsq1D, par=c(17.6, log10(1.7), log10(3), 17.4, log10(13), log10(0.7)),
rad=ellipses_box$ellipses$radhi[4:19], SB=ellipses_box$ellipses$SB[4:19], pixscale=0.2,
method='L-BFGS-B', lower=lower, upper=upper)$par

magplot(ellipses_box$ellipses$radhi[4:19], ellipses_nobox$ellipses$SB[4:19],
ylim=c(25,17), grid=TRUE, type='l')
points(ellipses_box$ellipses$radhi[4:19],ellipses_box$ellipses$SB[4:19])
#A simple bulge+disk surface brightness profile:
rlocs=seq(1,30,by=0.1)
bulge=profitRadialSersic(rlocs, mag=fit1D[1], re=10^fit1D[2], nser=10^fit1D[3])
disk=profitRadialSersic(rlocs, mag=fit1D[4], re=10^fit1D[5], nser=10^fit1D[6])
lines(rlocs, profoundFlux2SB(bulge, pixscale=0.2), col='red')
lines(rlocs, profoundFlux2SB(disk, pixscale=0.2), col='blue')
lines(rlocs, profoundFlux2SB(bulge+disk, pixscale=0.2), col='green')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ ellipse }% use one of  RShowDoc("KEYWORDS")
