\name{profoundInterp2d}
\alias{profoundInterp2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
2D image interpolation
}
\description{
A low level routine to interpolate a 2D image matrix at an arbitrary x/y pixel location. This function is unlikely to be used by the user, but it used internally to ensure that point sources defined by empirical point spread functions (PSFs) are accurately generated on an image.
}
\usage{
profoundInterp2d(x, y, image)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The x position at which to make the interpolation with respect to the x centre of \option{image}.
}
  \item{y}{
The x position at which to make the interpolation with respect to the x centre of \option{image}.
}
  \item{image}{
The image matrix to be used for the interpolation.
}
}
\details{
In practice this is a low level routine unlikely to be used by the user.

For this function (and really, it is for user ease when interpolating a PSF) [0,0] is always the R image centre of the input \option{image}. This means it would be at the usual [1.5,2] position of a 3x4 image matrix.
}
\value{
Matrix; a three column matrix where column 1 is the requested x interpolation locations, column 2 is the requested y interpolation locations and column 3 is the interpolated values.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundMakeSkyMap}}, \code{\link{profoundMakeSkyGrid}}
}
\examples{
# Nothing here
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ interpolate }
\keyword{ interpolation }
