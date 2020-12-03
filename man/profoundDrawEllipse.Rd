\name{profoundDrawEllipse}
\alias{profoundDrawEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Ellipse
}
\description{
Draws multiple ellipses on a plot window.
}
\usage{
profoundDrawEllipse(xcen = 0, ycen = 0, rad = 1, axrat = 1, ang = 0, box = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Numeric vector; x centre/s of the ellipse/s.
}
  \item{ycen}{
Numeric vector; y centre/s of the ellipse/s.
}
  \item{rad}{
Numeric vector; the major axis extent of the ellipse/s.
}
  \item{axrat}{
Numeric vector; the axial ratio of the ellipse/s as given by \option{radlo}/\option{radhi}.
}
  \item{ang}{
Numeric vector; the angle of the ellipse/s in the usual ProFit sense, see \code{profitMakeModel}.
}
  \item{box}{
Numeric vector; the boxiness of the ellipse/s in the usual ProFit sense, see \code{profitMakeModel}.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{lines}} to draw the ellipse/s.
}
}
\details{
This function uses all the standard \code{ProFit} conventions to define the input parameters
}
\value{
No value is returned, this function is run purely for the side effect of drawing an ellipse.

}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundGetEllipsesPlot}}, \code{\link{profoundGetEllipses}}, \code{\link{profoundGetEllipse}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, magzero=30, verbose=TRUE, plot=TRUE)
profoundDrawEllipse(profound$segstats$xcen, profound$segstats$ycen,
profound$segstats$R100/0.339, profound$segstats$axrat, profound$segstats$ang,
col='white', lty=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ ellipse }% use one of  RShowDoc("KEYWORDS")
