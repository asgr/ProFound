\name{profoundGetEllipsesPlot}
\alias{profoundGetEllipsesPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create diagnostic plot of estimated iso-photal ellipses
}
\description{
Generates a useful plot merging a rapidly changing colour mapping with the estimated ellipses.
}
\usage{
profoundGetEllipsesPlot(image = NULL, ellipses = NULL, segim = NULL, segID = 1,
segellipseID = "all", pixscale = 1, col = rep(rainbow(10, s = 0.5), 4), border = "auto",
lty = 'auto', lwd = 'auto', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{ellipses}{
Data.frame; the ellipse information, but in practice the \option{ellipse} list output of \code{\link{profoundGetEllipses}}.
}
  \item{segim}{
Integer matrix; optional, the segmentation map of the image. This matrix *must* be the same dimensions as \option{image}.
}
  \item{segID}{
Integer scalar; optional, the desired \option{segim} segment to extract from the \option{image}.
}
  \item{segellipseID}{
Integer vector; the segellipseID to be plotted. The default of 'all' will display all ellipses.
}
  \item{pixscale}{
Numeric scalar; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). This should only be used if the radii columns in \option{ellipses} have already been scaled by the pixel scale.
}
  \item{col}{
The colour palette to be used for the background \option{image}. The default is chosen to be high contrast, to make it easier to compare the computed ellipses with the underlying isophotes.
}
  \item{border}{
The colour of the ellipse border drawn by \code{draw.ellipse}. If 'auto' then a sensible default is chosen.
}
  \item{lty}{
The line type of the ellipse border drawn by \code{draw.ellipse}. If 'auto' then a sensible default is chosen (\option{lty}=1 within the 90\% flux radius and \option{lty}=2 outside).
}
  \item{lwd}{
The line width of the ellipse border drawn by \code{draw.ellipse}. If 'auto' then a sensible default is chosen (\option{lwd}=0.5 within the 50\% flux radius, \option{lwd}=1 above the 50\% flux radius, except for the annuli at 50\%/90\% which is \option{lwd}=2).
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}.
}
}
\details{
The default options should create useful diagnostics, but there are lots of potential plots that can be made with the outputs of \code{\link{profoundGetEllipses}}, including e.g. making plots of how various parameters behave with radius, which can give helpful insight to starting parameters for bulge and disk profiles. The user is encouraged to experiment.
}
\value{
No value is returned, this function is run purely for the side effect of making a diagnostic plot.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundGetEllipses}}, \code{\link{profoundGetEllipse}}, \code{\link{profoundDrawEllipse}}
}
\examples{
\dontrun{
# We need the ProFit library to show the profile: library(ProFit)
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFit"))$imDat
segim = readFITS(system.file("extdata", 'KiDS/G266035segim.fits', package="ProFit"))$imDat
ellipses = profoundGetEllipses(image=image, segim=segim, segID=4, plot=FALSE)

#We can get a good overall idea of how good the ellipses are by running with defaults:
profoundGetEllipsesPlot(image=image, ellipses=ellipses$ellipses)

#We can check a specific ellipse too:
profoundGetEllipsesPlot(image=ellipses$segellipses==8, ellipses=ellipses$ellipses,
segellipseID=8, col=grey(0:1), border='red', lwd=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ ellipse }% use one of  RShowDoc("KEYWORDS")
