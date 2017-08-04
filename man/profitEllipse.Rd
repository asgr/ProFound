\name{profitEllipse}
\alias{profitEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure Isophotal Flux for Pseudo-Ellipses
}
\description{
In the world of galaxy fitting, projected 1D flux intensity (or surface brightness) plots are popular. This function implements the low level functionality of deprojecting image pixels given a set of geometrical parameters. In simple terms this means ellipses are expanded back up to circles. We use the term pseudo-ellipse since we can also account for boxiness distortion (if desired).
}
\usage{
profitEllipse(x, y, flux, xcen = 0, ycen = 0, ang = 0, axrat = 1, box = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Either a vector of image pixel midpoints (given in the usual ProFit standard, where pixel mid-points are half-integer), or an image matrix, which is then used for \option{x}, \option{y} and \option{flux}.
}
  \item{y}{
Vector; image pixel midpoints (given in the usual ProFit standard, where pixel mid-points are half-integer). Not required if an image matrix is being parsed to \option{x}.
}
  \item{flux}{
Vector; image pixel fluxes. Not required if an image matrix is being parsed to \option{x}.
}
  \item{xcen}{
Scalar; x centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{ycen}{
Scalar; y centre of the 2D Sersic profile (can be fractional pixel positions).
}
  \item{ang}{
Scalar; the orientation of the major axis of the Sersic profile in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang.
}
  \item{axrat}{
Scalar; axial ratio of the Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{box}{
Scalar; the boxiness of the Sersic profile that traces contours of iso-flux, defined such that r[mod]=(x^(2+box)+y^(2+box))^(1/(2+box)). When box=0 the iso-flux contours will be normal ellipses, but modifications between -1<box<1 will produce visually boxy distortions. Negative values have a pin-cushion effect, whereas positive values have a barrel effect (the major and minor axes staying fixed in all cases).
}
}
\details{
This function mostly exists for usage by the higher level \code{\link{profitEllipsePlot}} function. However, the outputs might well of interest if users want to create their own bespoke 1D plotting outputs.
}
\value{
A two column matrix containing columns rad (the effective projected major axis radius for the pixel) and flux (the un-corrected flux for this pixel).
}
\references{
As noted above, the pixel flux is not corrected for distortion. This is rarely an issue, but for highly elliptical galaxies where the gradient changes radically through the pixel along the minor versus major axis the isophotal contour will be biased. In practice this bias will be the same as for the model, so pixel to pixel comparisons are still fairly valid. If this level of inconsistency annoys you, then you almost certainly should not be attempting to make deprojected 1D plots of PSF convolved bulge-disk systems. This is because once a pure elliptical disk has been convolved with a PSF it is not strictly possible to project this to circular annuli (i.e. the pseudo major-axis output we desire) using elliptical isophotes. C'est la vie.

By ProFit convention the bottom-left part of the bottom-left pixel when plotting the image matrix is c(0,0) and the top-right part of the bottom-left pixel is c(1,1), i.e. the mid-point of pixels are half integer values in x and y.

To confuse things a bit, when R plots an image of a matrix it is transposed and re-ordered vertically to how it appears if you print the matrix directly to screen, i.e. compare print(matrix(1:4,2,2)) and image(matrix(1:4,2,2)). The lowest value (1) is top-left when printed but bottom-left when displayed using image (the red pixel). Both are "correct": the issue is whether you consider the first element of a matrix to be the Cartesian x position (movement in x) or a row element (movement in y). Matrices in maths are always written top-left first where the first argument refers to row number, but images by convention are accessed in a Cartesian sense. Hence [3,4] in a maths matrix means 3 down and 4 right from the top-left, but 3 right and 4 up from the bottom-left in an image.
}
\author{
Aaron Robotham
}
\note{
Projecting data back to 1D profiles is a knowingly imperfect process, but it can be useful when exploring the data and for producing explanatory plots. In most use cases the data error and systematics will probably dominate the process, but be cautious in overly relying on these outputs- the 2D fit residuals from \code{\link{profitMakePlots}} contains more information. Where the 1D and 2D plots appear to be in conflict re the quality of the fit, the latter should always take precedence.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profitEllipsePlot}}, \code{\link{profitRemakeModellist}}
}
\examples{

#The rough best fit model for G266033 (KiDS)

model=list(
  sersic=list(
    xcen = c(65.60642, 65.60642),
    ycen = c(78.6091, 78.6091),
    mag = c(18.49816, 16.97754),
    re = c(1.69112, 14.75940),
    nser = c(1.053961, 1),
    ang = c(39.53360, 35.02479),
    axrat = c(1, 0.5990869),
    box = c(0,0)
  )
)

data('ExampleInit')
image = readFITS(system.file("extdata", 'KiDS/G266033fitim.fits',package="ProFound"))$imDat

temp=profitEllipse(x=image, xcen=model$sersic$xcen[2], ycen=model$sersic$ycen[2], ang=
  model$sersic$ang[2], axrat=model$sersic$axrat[2], box=model$sersic$box[2])

#A rough deprojected ellipse plot:

magplot(temp, type='l', log='y', xlim=c(0,50), ylim=c(1e-12,2e-9), xlab='Pixels',
 ylab='Pixel Flux')

#Notice in the example that the core looks noisier, this is because we are deprojecting
#the bulge using the disk parameters. We can do the same using the bulge:

temp2=profitEllipse(x=image, xcen=model$sersic$xcen[1], ycen=model$sersic$ycen[1], ang=
  model$sersic$ang[1], axrat=model$sersic$axrat[1], box=model$sersic$box[1])
  
magplot(temp2, type='l', log='y', xlim=c(0,50), ylim=c(1e-12,2e-9), xlab='Pixels',
 ylab='Pixel Flux')

#The inner region (<5 pixels) is better deprojected using the bulge parameters.

#There is no simple way to meaningfully deproject such systems, sincew there will always
#be regions that are an even mix of bulge and disk flux, and these require different
#deprojection parameters. This is even true with IRAF Ellipse style ring fitting. The
#collapsing of data in this manner is an inherently lossy process!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ellipse }
