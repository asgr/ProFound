\name{profoundCover}
\alias{profoundCover}
\alias{profoundAperCover}
\alias{profoundEllipCover}
\alias{profoundPolyCover}
\alias{profoundAperFlux}
\alias{profoundEllipFlux}
\alias{profoundPolyFlux}
\alias{profoundSegimFlux}
\alias{profoundAperWeight}
\alias{profoundEllipWeight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fraction of Pixels in Shapes
}
\description{
Calculates the fraction of pixels inside an aperture, ellipse or polygon using recursive descent in Rcpp.
}
\usage{
profoundAperCover(x, y, cx, cy, rad, depth = 3L, nthreads = 1L)
profoundEllipCover(x, y, cx, cy, rad, ang = 0, axrat = 1, depth = 3L, nthreads = 1L)
profoundPolyCover(x, y, poly_x, poly_y, depth = 3L, nthreads = 1L)

profoundAperFlux(image, cx, cy, rad, wt = as.numeric(c(1)), rad_re = as.numeric(c(1),
  deblend = FALSE, depth = 3L, nthreads = 1L)
profoundEllipFlux(image, cx, cy, rad, ang = as.numeric(c(0)), axrat = as.numeric(c(1)),
  wt = as.numeric(c(1)), rad_re = as.numeric(c(1), deblend = FALSE, depth = 3L,
  nthreads = 1L)
profoundPolyFlux(image, poly_x, poly_y, depth = 3L, nthreads = 1L)
profoundSegimFlux(image, segim, nthreads = 1L)

profoundAperWeight(cx, cy, rad, dimx = 100L, dimy = 100L, wt = as.numeric(c(1)),
  rad_re = as.numeric(c(1), depth = 3, nthreads = 1)
profoundEllipWeight(cx, cy, rad, ang = as.numeric(c(0)), axrat = as.numeric(c(1)),
  dimx = 100L, dimy = 100L, wt = as.numeric(c(1)), rad_re = as.numeric(c(1), depth = 3,
  nthreads = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the pixel x centre positions (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{y}{
Numeric vector; the pixel y centre positions (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{image}{
Numeric matrix; the image we want to compute flux on.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{mask}.
}
  \item{cx}{
Numeric scalar/vector; the x centre of the aperture (in pixels, where pixel centres are half integer in \code{ProFound}). Scalar for 'Cover' functions and vector for 'Flux' and 'Weight' functions.
}
  \item{cy}{
Numeric scalar/vector; the y centre of the aperture (in pixels, where pixel centres are half integer in \code{ProFound}). Scalar for 'Cover' functions and vector for 'Flux' and 'Weight' functions (must be same length as \option{cx}).
}
  \item{rad}{
Numeric scalar/vector; the radius or semi-major axis of the aperture / ellipse (in pixels). Scalar for 'Cover' functions and vector for 'Flux' and 'Weight' functions (in these cases, if length 1 then value is repeated for all entries, otherwise must be same length as \option{cx}).
}
  \item{ang}{
Numeric scalar; orientation of the semi-major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Scalar for 'Cover' functions and vector for 'Flux' and 'Weight' functions (in these cases, if length 1 then value is repeated for all entries, otherwise must be same length as \option{cx}).
}
  \item{axrat}{
Numeric scalar; the semi-minor/semi-major axial ratio of the ellipse (bounded between 0 - 1). Scalar for 'Cover' functions and vector for 'Flux' and 'Weight' functions (in these cases, if length 1 then value is repeated for all entries, otherwise must be same length as \option{cx}).
}
  \item{wt}{
Numeric scalar/vector; the per aperture / ellipse weight to apply when computing the weight map and/or fluxes. If length 1 then value is repeated for all entries, otherwise must be same length as \option{cx}.
}
  \item{rad_re}{
Numeric scalar; the exponential Re scaling to apply to the weights. Setting to 0 means the weight response is flat throughout the aperture / ellipse. Setting to the true effective radius of the profile will create a more accurate deblend, but the minimum value should reflect the PSF FWHM (where Re = FWHM/2).
}
  \item{dimx}{
Integer scalar; the number of pixels in the x-dimension for the resulting weight map.
}
  \item{dimy}{
Integer scalar; the number of pixels in the y-dimension for the resulting weight map.
}
  \item{poly_x}{
Numeric vector; x vertices of polygon. This need not be closing, so 4 points define square not 5.
}
  \item{poly_y}{
Numeric vector; y vertices of polygon (must be same length as \option{poly_x}). This need not be closing, so 4 points define square, not 5.  
}
  \item{deblend}{
Logical; should an internal weight map (made via the associated Weight functions) be used to normalise the flux across competing apertures/ellipses. This allows for basic deblending of flux.  
}
  \item{depth}{
Integer scalar; the recursion depth to use. Depth 3 (default) is a good trade off of accuracy and speed (see Examples), with 4 probably being as high as you would want to set it for a typical ues-case.
}
  \item{nthreads}{
Integer scalar; the number of CPU threads to run on. Only valid if \code{ProFound} is built with OpenMP support.
}
}
\details{
Re the recursion \option{depth}, higher numbers are more accurate, where you probably do not want to set it below 3 (the default). Above 5 will suffer notable slow-down. 4 is usually about the performance versus accuracy sweet spot. Note SExtractor uses simple over-sampling instead, where the default is a factor 5. Since 2^3 = 8 and 2^2 = 4, this is more like \option{depth} = 2 (2.32 to be precise).

When \option{rad_re} is not 0 an exponential light profile is assumed for deblending, where the value of \option{rad_re} is the half light radius (AKA effective radius) in pixels. This works exactly when the object in question has an exponential light profile (a disk galaxy, say) but also works reasonably for PSF sources and steeper profile objects (elliptical galaxies).

When \option{rad_re} is not 0 and a full vector of \option{wt} is supplied, the specific \option{wt} values are divided by \option{rad_re}^2 (aperture) or \option{rad_re}^2 x \option{axrat} (ellipse). This is reflect the idea that the \option{wt} is the total flux expected within the aperture, so a different area of aperture/ellipse will imply a different surface brightness. Scaling like this means we get better flux deblending when we have large differences in aperture size.

\code{profoundAperFlux} and \option{profoundAperFlux} can be run iteratively, where the flux output of a run becomes the \option{wt} input for a new run. Generally running it like this will create converged flux in 3-5 steps.

These functions are written in performant Rcpp, and can be built against OpenMP where available (setting \option{nthreads} higher than 1 if appropriate).
}
\value{
For the Cover functions: Numeric vector the same length as \option{x} and \option{y}, which is what fraction of the pixel lies under the desired aperture / ellipse.

For the Flux functions: Numeric scalar, the computed image flux within the region requested. This uses the "cover" fractional pixel weights when \option{depth} is larger than 0.

\code{profoundSegimFlux} is a special case, wehre the output is a vector the same length as the maximum segID on \option{segim}, where each value corresponds to the summed flux of that segment by index, i.e. output[1] is the total flux in segment 1, and output[N] is the total flux in segment N etc. Some might be 0 in cases where that segID is not present at all in the \option{segim}.

For the Weight functions. Numeric matrix; the total number of contributing apertures/ellipses. This can be non-integer when \option{depth} is above 0. Where valued 1 all of a pixels flux belongs to a single aperture, and value 2 it would be shared between two (50\% each) etc. This is mostly used as an optional input the Flux functions (which can take \option{weight} as in input).
}
\references{
\code{profoundPolyCover} uses code provided on \url{https://wrfranklin.org/Research/Short_Notes/pnpoly.html} by W. Randolph Franklin.
}
\author{
Aaron Robotham and Rodrigo Tobar
}
\seealso{
\code{\link{profoundAperPhot}}, \code{\link{profoundEllipseSeg}}
}
\examples{
library(magicaxis)
library(plotrix)

temp = expand.grid(1:10 - 0.5, 1:10 - 0.5)

output = matrix(profoundAperCover(temp[,1], temp[,2], 5, 5, 3, depth=4), 10, 10)
magimage(output, magmap=FALSE, zlim=c(0,1))
draw.circle(5, 5, radius=3.2, border='red')

output = matrix(profoundAperCover(temp[,1], temp[,2], 3.2, 3.2, 1.8, depth=4), 10, 10)
magimage(output, magmap=FALSE, zlim=c(0,1))
draw.circle(3.2, 3.2, radius=1.8, border='red')

output = matrix(profoundAperCover(temp[,1], temp[,2], 7.5, 7.5, 5.8, depth=4), 10, 10)
magimage(output, magmap=FALSE, zlim=c(0,1))
draw.circle(7.5, 7.5, radius=5.8, border='red')

#Testing accuracy:

\dontrun{
ap_rad = 0.5 #this is smaller than a single pixel, so hard test case
ap_area = pi*ap_rad^2

ap_grid = expand.grid(seq(0,1,by=0.01), seq(0,1,by=0.01))
pix_grid = expand.grid(seq(-0.5,1.5,by=1), seq(-0.5,1.5,by=1))

test = {}
for(i in 1:dim(ap_grid)[1]){
  test = c(test, sum(profoundAperCover(pix_grid[,1], pix_grid[,2], ap_grid[i,1],
  ap_grid[i,2], ap_rad, depth=3))/ap_area)
}

#with depth 3 we do not do worse than 6.5\% accuracy:
magimage(matrix(test,101,101), magmap=FALSE, zlim=c(0.9,1.1))
maghist(test)
}

#Some ellipse examples:

output = matrix(profoundEllipCover(temp[,1], temp[,2], 5, 5, 4.5, 30, 0.5, depth=4),
  10, 10)
magimage(output, magmap=FALSE, zlim=c(0,1))
draw.ellipse(5, 5, 4.5, 4.5*0.5, 30+90, border='red') # +90 because different ang def

output = matrix(profoundEllipCover(temp[,1], temp[,2], 4.2, 5.1, 2.5, 110, 0.8, depth=4),
  10, 10)
magimage(output, magmap=FALSE, zlim=c(0,1))
draw.ellipse(4.2, 5.1, 2.5, 2.5*0.8, 110+90, border='red') # +90 because different ang def

#poly example:

poly_x = c(0.958, 3.296, 1.583, 5.940, 8.969, 9.886, 5.701, 3.575, 2.964)
poly_y = c(1.482, 4.677, 8.059, 4.977, 9.149, 2.647, 0.016, 2.798, 2.535)

output = matrix(profoundPolyCover(temp[,1], temp[,2], poly_x, poly_y, depth=4),10,10)
magimage(output, magmap=FALSE, zlim=c(0,1))
polygon(temp,border='red')

#Weight map:
set.seed(666)
ran_x = runif(20,0,100)
ran_y = runif(20,0,100)
ran_rad = runif(20,3,20)

magimage(profoundAperWeight(ran_x, ran_y, ran_rad), magmap=FALSE)
for(i in seq_along(ran_x)){
  draw.circle(ran_x[i], ran_y[i], ran_rad[i],, border='red')
}

ran_ang = runif(20,0,180)
ran_axrat = runif(20,0,1)

magimage(profoundEllipWeight(ran_x, ran_y, ran_rad, ran_ang, ran_axrat), magmap=FALSE)
for(i in seq_along(ran_x)){
  draw.ellipse(ran_x[i], ran_y[i], ran_rad[i], ran_rad[i]*ran_axrat[i], ran_ang[i]+90,
  border='red')
}

#Provide rad_re. This will bias flux towards smaller isophotal contour sources.
magimage(profoundEllipWeight(ran_x, ran_y, ran_rad, ran_ang, ran_axrat, rad_re=10),
  magmap=TRUE)
for(i in seq_along(ran_x)){
  draw.ellipse(ran_x[i], ran_y[i], ran_rad[i], ran_rad[i]*ran_axrat[i], ran_ang[i]+90,
  border='red')
}

\dontrun{
#Some elliptical photometry:

image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
                                     package="ProFound"), header=FALSE)

profound = profoundProFound(image, verbose=TRUE, plot=TRUE)
ellip_ref = profound$segstats[,c("flux", "xcen", "ycen", "R50", "R100",  "ang", "axrat")]

#with rad_re set you can see the deblend solution that will be used for splitting flux:
magimage(profoundEllipWeight(ellip_ref$xcen, ellip_ref$ycen, ellip_ref$R100,
  ellip_ref$ang, ellip_ref$axrat, 356, 356, rad_re=ellip_ref$R50), magmap=TRUE, add=TRUE,
  col=c(NA, hsv(1, s=(0:1000)/1000, alpha=0.5)))
  
ellip_flux_noblend =  profoundEllipFlux(image, ellip_ref$xcen, ellip_ref$ycen,
  ellip_ref$R100, ellip_ref$ang, ellip_ref$axrat, deblend=FALSE)
ellip_flux_deblend =  profoundEllipFlux(image, ellip_ref$xcen, ellip_ref$ycen,
  ellip_ref$R100, ellip_ref$ang, ellip_ref$axrat, deblend=TRUE, rad_re=ellip_ref$R50)
  
magplot(profound$segstats$flux, profound$segstats$flux/ellip_flux_noblend, log='xy',
  xlab='ProFound [Flux]', ylab='ProFound / Ellip [Flux]', ylim=c(0.1,10))
points(profound$segstats$flux, profound$segstats$flux/ellip_flux_deblend, col='red')
legend('topleft', legend=c('No Deblend', 'Deblend'),
  col=c('black', 'red'), pch=1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
