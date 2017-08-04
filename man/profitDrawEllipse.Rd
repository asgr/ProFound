\name{profitDrawEllipse}
\alias{profitDrawEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Ellipse
}
\description{
Calculates the \option{x} and \option{y} location of an ellipse, allowing for the presence of boxiness.
}
\usage{
profitDrawEllipse(xcen = 0, ycen = 0, rad = 1, axrat = 1, ang = 0, box = 0, ...)
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
Numeric vector; the angle of the ellipse/s in the usual ProFit sense, see \code{\link{profitMakeModel}}.
}
  \item{box}{
Numeric vector; the boxiness of the ellipse/s in the usual ProFit sense, see \code{\link{profitMakeModel}}.
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
\code{\link{profitGetEllipsesPlot}}, \code{\link{profitGetEllipses}}, \code{\link{profitGetEllipse}}
}
\examples{
\dontrun{
library(ProFit)

model = list(
 sersic = list(
  xcen = 50,
  ycen = 50,
  mag = 15,
  re = 10,
  nser = 4,
  ang = 30,
  axrat = 0.3,
  box = 0.5
 )
)

image=profitMakeModel(model)$z
temp=magimage(image)
contour(temp, col='red', drawlabels=FALSE, add=TRUE)
profitDrawEllipse(model$sersic$xcen, model$sersic$ycen, rad=26, axrat=model$sersic$axrat,
ang=model$sersic$ang, box=model$sersic$box, col='blue', lty=2)

model = list(
 sersic = list(
  xcen = 50,
  ycen = 50,
  mag = 15,
  re = 10,
  nser = 4,
  ang = 30,
  axrat = 0.3,
  box = -0.5
 )
)

image=profitMakeModel(model)$z
temp=magimage(image)
contour(temp, col='red', drawlabels=FALSE, add=TRUE)
profitDrawEllipse(model$sersic$xcen, model$sersic$ycen, rad=30, axrat=model$sersic$axrat,
ang=model$sersic$ang, box=model$sersic$box, col='blue', lty=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ellipse }% use one of  RShowDoc("KEYWORDS")
